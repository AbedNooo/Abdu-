import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import truncnorm
from scipy.special import expit

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
sns.set_style("whitegrid")
plt.rcParams.update({'font.size': 12, 'figure.titlesize': 16, 'figure.titleweight': 'bold'})

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================
def draw_truncated_normal(mu, sigma, lower, upper, size, rng):
    a, b = (lower - mu) / sigma, (upper - mu) / sigma
    return truncnorm.rvs(a, b, loc=mu, scale=sigma, size=size, random_state=rng)

# ==============================================================================
# 3. CLASS DEFINITIONS
# ==============================================================================
class SolidarityFund:
    def __init__(self, initial_capital, reserve_target=25000.0): 
        self.balance = float(initial_capital)
        self.Q_outstanding = 20000.0          # Qard Hasan (Debt)
        self.reserve_target = float(reserve_target)

    def process_year(self, subsidies, weather_shocks, expected_losses, current_tau, adoption_mask, trigger=0.70):
        # 1. INFLOW
        contributions_inflow = np.sum(subsidies[adoption_mask] * float(current_tau))
        self.balance += contributions_inflow

        # 2. OUTFLOW (Pro-Rata Fairness)
        # Trigger is now based on deviations from Normal (1.0). e.g., < 0.75
        triggered_indices = np.where((weather_shocks < trigger) & (adoption_mask))[0]
        payouts = np.zeros(len(subsidies), dtype=float)
        
        if triggered_indices.size > 0:
            total_claimed_loss = np.sum(expected_losses[triggered_indices])
            
            if self.balance >= total_claimed_loss:
                payouts[triggered_indices] = expected_losses[triggered_indices]
                self.balance -= total_claimed_loss
            else:
                coverage_ratio = self.balance / total_claimed_loss
                payouts[triggered_indices] = expected_losses[triggered_indices] * coverage_ratio
                self.balance = 0.0 

        # 3. EARLY REPAYMENT RULE
        safety_floor = 15000.0 
        if self.Q_outstanding > 0 and self.balance > safety_floor:
            early_pay = min(0.10 * contributions_inflow, self.Q_outstanding)
            self.Q_outstanding -= early_pay
            self.balance -= early_pay

        # 4. SURPLUS LOGIC (Debt First)
        surplus = max(0.0, self.balance - self.reserve_target)
        dividend_per_agent = 0.0

        if surplus > 0:
            if self.Q_outstanding > 0:
                repayment = min(surplus, self.Q_outstanding)
                self.Q_outstanding -= repayment
                self.balance -= repayment
                surplus -= repayment  

            if self.Q_outstanding <= 1.0 and surplus > 0:
                num_members = int(np.sum(adoption_mask))
                if num_members > 0:
                    distributable = surplus * 0.5 
                    dividend_per_agent = distributable / num_members
                    self.balance -= distributable

        return payouts, dividend_per_agent

class Agent:
    def __init__(self, agent_id, params, rng):
        self.id = agent_id
        # Small Farm Reality
        raw_area = rng.lognormal(mean=1.42, sigma=0.5)
        self.Area = float(np.clip(raw_area, 0.5, 50.0))

        # Liquidity Trap
        self.W = float(rng.uniform(300, 600) * self.Area) 
        
        # Subsidy: Conditional Eco-Scheme Payment
        self.subsidy_amt = float(174.0 * self.Area)  
        
        # Recalibrated Yield Potential:
        # Since avg Omega is now 1.0 (was 0.7), we lower Y_max to keep revenue realistic.
        self.Y_max = float(4000.0 * self.Area)       
        
        self.rho = float(params['rho'])     
        self.beta0 = float(params['beta0']) 
        
        self.adopted = 0
        self.t_adopt = -1 
        
        self.N_curr = 150.0 
        self.a = 0.05
        self.c_n = 0.80; self.chi = 0.005
        self.a0 = 15.0; self.c_op = 10.0; self.c_comm = 5.0; self.c_ha = 5.0

        # History
        self.nue_history = []
        self.target_history = []

    def calculate_mrv_cost(self):
        return (self.a0 + self.c_op + self.c_comm) + (self.c_ha * self.Area)

# ==============================================================================
# 4. SINGLE SIMULATION RUN
# ==============================================================================
def run_single_simulation(seed, use_fund=True):
    rng = np.random.default_rng(seed)
    n_agents = 200
    t_max = 15

    # --- CALIBRATION ---
    beta1 = 0.015  # Money
    
    # Institutional Cohesion (Peer Effect)
    beta2_base = 2.0  
    cohesion_k = 0.35 
    if use_fund:
        beta2 = beta2_base * (1.0 + cohesion_k) 
    else:
        beta2 = beta2_base                      
    
    beta3 = -2.5   # Fear

    # --- WEATHER RE-CALIBRATION ---
    # Normal Year = 1.0 (100% Potential).
    # Drought < 0.75.
    # This allows Target = 0.60 * Omega to yield ~60% average.
    omega_mu = 1.0
    omega_sigma = 0.15
    omega_lower = 0.50  # Severe Drought
    omega_upper = 1.30  # Bumper Crop
    
    # Trigger Payout if yield drops below 75% of potential (25% loss)
    trigger_omega = 0.75 

    agents = []
    beta0_pool = np.concatenate([rng.normal(-5.0, 1.0, int(n_agents*0.74)), rng.normal(-20.0, 1.0, int(n_agents*0.26))])
    rng.shuffle(beta0_pool)

    for i in range(n_agents):
        p = {'rho': float(np.clip(rng.normal(0.6, 0.1), 0.1, 0.9)), 'beta0': float(beta0_pool[i])}
        agents.append(Agent(i, p, rng))

    fund = SolidarityFund(initial_capital=20000.0, reserve_target=25000.0)
    
    adoption_series = []
    fund_series = []
    debt_series = []
    nue_series = []

    for t in range(t_max):
        # --- FIX 2: SYSTEMIC REGIONAL WEATHER ---
        # One main shock for the region
        omega_t = draw_truncated_normal(mu=omega_mu, sigma=omega_sigma, lower=omega_lower, upper=omega_upper, size=1, rng=rng)[0]
        # Small local variations (micro-climates)
        omega = np.clip(omega_t + rng.normal(0.0, 0.05, n_agents), omega_lower, omega_upper)
        
        current_tau = 0.25 if t < 3 else 0.15
        prev_adoption_rate = adoption_series[-1] / 100.0 if t > 0 else 0.0
        
        subsidies = np.array([a.subsidy_amt for a in agents])
        adoption_mask = np.array([a.adopted == 1 for a in agents], dtype=bool)
        areas = np.array([a.Area for a in agents])
        exp_losses = 30.0 * areas 
        
        # Fund Processing
        if use_fund:
            payouts, dividend = fund.process_year(subsidies, omega, exp_losses, current_tau, adoption_mask, trigger_omega)
        else:
            payouts = np.zeros(n_agents); dividend = 0.0

        current_adoption = 0
        sys_nue_num = 0.0; sys_nue_den = 0.0

        for i, ag in enumerate(agents):
            # Agronomy
            N_applied = ag.N_curr
            Y_total = ag.Y_max * (1.0 - np.exp(-ag.a * N_applied)) * omega[i]
            N_uptake = 0.018 * Y_total
            
            # Conditionality
            true_nue = (N_uptake / (N_applied * ag.Area)) if N_applied > 0 else 0.0
            obs_nue = np.clip(true_nue + rng.normal(0.0, 0.03), 0.0, 1.0)
            
            # --- FIX 1: THESIS-CONSISTENT TARGET ---
            # Target is 60% of Realized Potential.
            # Since avg Omega ~ 1.0, avg Target ~ 0.60 (60%).
            target_val = 0.60 * omega[i]

            ag.nue_history.append(obs_nue)
            ag.target_history.append(target_val)
            if len(ag.nue_history) > 3: ag.nue_history.pop(0)
            if len(ag.target_history) > 3: ag.target_history.pop(0)

            clawback = 0.0
            # Grace Period (2 years from INDIVIDUAL adoption)
            if ag.adopted and ag.t_adopt != -1 and (t - ag.t_adopt) >= 2:
                avg_perf = float(np.mean(ag.nue_history))
                avg_target = float(np.mean(ag.target_history))
                if avg_perf < avg_target:
                    shortfall = avg_target - avg_perf
                    # Cap Fix: Limited by Subsidy amount
                    clawback = min(shortfall * 3.0 * ag.subsidy_amt, ag.subsidy_amt)

            # Economics
            cost_N = ((ag.c_n * N_applied) + (0.5 * ag.chi * (N_applied**2))) * ag.Area
            residual_ops_cost = 600.0 * ag.Area
            cost_mrv_actual = ag.calculate_mrv_cost() if ag.adopted else 0.0
            revenue = 0.35 * Y_total 
            
            net_subsidy_cash = 0.0
            if ag.adopted:
                sol_contrib = (current_tau * ag.subsidy_amt) if use_fund else 0.0
                net_subsidy_cash = ag.subsidy_amt - sol_contrib - cost_mrv_actual - clawback

            profit = (revenue - cost_N - residual_ops_cost) + net_subsidy_cash + payouts[i] + (dividend if (ag.adopted and use_fund) else 0)
            
            # Wealth Update
            if profit > 0: ag.W += profit * 0.20 
            else: ag.W += profit 

            # ADOPTION LOGIC
            potential_mrv_total = ag.calculate_mrv_cost()
            potential_mrv_per_ha = potential_mrv_total / ag.Area
            
            upfront_cost = 80.0 + (potential_mrv_total * 0.5)
            feasibility = 1.0 if ag.W >= upfront_cost else 0.0

            exp_contrib_per_ha = (current_tau * 174.0) if use_fund else 0.0
            net_incentive_per_ha = 174.0 - exp_contrib_per_ha - potential_mrv_per_ha
            
            trust_mult = 1.0 
            coverage_ratio = ag.W / (upfront_cost + 1.0)
            trust_mult = 1.0 + (1.5 - np.clip(coverage_ratio, 0.5, 2.5)) * 0.2
            
            rho_effective = ag.rho * trust_mult
            if use_fund: rho_effective *= 0.20

            exponent = ag.beta0 + beta1 * net_incentive_per_ha + beta2 * prev_adoption_rate + beta3 * rho_effective
            p_adopt = expit(exponent) * feasibility 

            if (not ag.adopted) and (rng.random() < p_adopt): 
                ag.adopted = 1
                ag.t_adopt = t 

            # Learning
            target_N = 105.0 if ag.adopted else 150.0
            ag.N_curr = (1.0 - 0.40) * ag.N_curr + 0.40 * target_N

            if ag.adopted: current_adoption += 1
            sys_nue_num += N_uptake; sys_nue_den += (N_applied * ag.Area)

        adoption_series.append(current_adoption / n_agents * 100.0)
        nue_series.append((sys_nue_num / sys_nue_den) * 100.0 if sys_nue_den > 0 else 0.0)
        fund_series.append(fund.balance)
        debt_series.append(fund.Q_outstanding)

    return adoption_series, fund_series, debt_series, nue_series

# ==============================================================================
# 5. MONTE CARLO WRAPPER
# ==============================================================================
def run_monte_carlo(n_runs=100, use_fund=True):
    all_adopt = []
    all_fund = []
    all_debt = []
    all_nue = []
    
    print(f"Running {n_runs} simulations (Fund={use_fund})...")
    for i in range(n_runs):
        a, f, d, n = run_single_simulation(seed=i, use_fund=use_fund)
        all_adopt.append(a)
        all_fund.append(f)
        all_debt.append(d)
        all_nue.append(n)
        
    return np.array(all_adopt), np.array(all_fund), np.array(all_debt), np.array(all_nue)

# ==============================================================================
# 6. VISUALIZATION
# ==============================================================================
if __name__ == "__main__":
    runs = 100 
    
    # Run
    adopt_F, fund_F, debt_F, nue_F = run_monte_carlo(runs, use_fund=True)
    adopt_NF, _, _, nue_NF = run_monte_carlo(runs, use_fund=False)

    # Plotter
    def plot_smooth(ax, data, color, label, linestyle='-'):
        mean = np.mean(data, axis=0)
        lower = np.percentile(data, 10, axis=0)
        upper = np.percentile(data, 90, axis=0)
        ax.plot(mean, color=color, linewidth=3, linestyle=linestyle, label=label)
        ax.fill_between(range(len(mean)), lower, upper, color=color, alpha=0.15)

    fig, axs = plt.subplots(1, 3, figsize=(18, 6))

    # A. Adoption
    plot_smooth(axs[0], adopt_F, '#2ca02c', 'With Fund')
    plot_smooth(axs[0], adopt_NF, 'gray', 'Without Fund', '--')
    axs[0].set_title(f'A. Adoption Rate (Mean of {runs} runs)', fontweight='bold')
    axs[0].set_ylabel('% Farmers')
    axs[0].legend(loc='upper left')
    axs[0].grid(True, alpha=0.3)

    # B. Fund Health
    plot_smooth(axs[1], fund_F, 'green', 'Fund Balance')
    plot_smooth(axs[1], debt_F, 'red', 'Qard Hasan (Debt)', '--')
    axs[1].axhline(25000, ls=':', color='gray', label='Reserve Target')
    axs[1].set_title('B. Fund Solvency & Repayment', fontweight='bold')
    axs[1].set_ylabel('EUR')
    axs[1].legend(loc='center right')
    axs[1].grid(True, alpha=0.3)

    # C. NUE
    plot_smooth(axs[2], nue_F, 'purple', 'With Fund')
    plot_smooth(axs[2], nue_NF, 'gray', 'Without Fund', '--')
    axs[2].set_title('C. Nitrogen Efficiency', fontweight='bold')
    axs[2].set_ylabel('NUE (%)')
    axs[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
