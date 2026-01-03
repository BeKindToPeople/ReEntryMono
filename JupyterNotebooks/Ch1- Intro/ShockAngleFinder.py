import math

def calculate_oblique_shock_beta(Mach, theta_deg, gamma=1.4):
    """
    Solves the Theta-Beta-Mach relation for beta (shock angle)
    using the Newton-Raphson method.
    """
    
    # 1. Convert everything to Radians immediately
    # Dimensions: [theta] = rad
    theta = math.radians(theta_deg)
    
    # Check for detached shock limit (simplified check)
    # If theta is too high, no solution exists. 
    # For M=36, theta=15 is safe.
    
    # 2. Define the Target Function f(beta) = 0
    # derived from: tan(theta) = 2 cot(beta) * ...
    def theta_beta_m_func(beta):
        # Prevent division by zero if beta is 0 or 90
        if beta <= 0 or beta >= math.pi/2:
            return 1e9 
            
        numerator = (Mach**2 * math.sin(beta)**2 - 1)
        denominator = (Mach**2 * (gamma + math.cos(2 * beta)) + 2)
        
        # Calculate the Right Hand Side (RHS)
        rhs = 2 * (1 / math.tan(beta)) * (numerator / denominator)
        
        # We want RHS - tan(theta) = 0
        return rhs - math.tan(theta)

    # 3. Initial Guess for Beta
    # Physical reasoning:
    # The Weak Shock solution is usually slightly larger than theta.
    # The Strong Shock solution is close to 90 degrees.
    # We want the Weak Shock.
    beta_guess = math.radians(theta_deg + 2.0) 
    
    # 4. Newton-Raphson Loop
    tolerance = 1e-6
    max_iterations = 20
    epsilon = 1e-6 # Step size for numerical derivative
    
    print(f"{'Iter':<5} | {'Beta (deg)':<12} | {'Error':<12}")
    print("-" * 35)

    for i in range(max_iterations):
        f_val = theta_beta_m_func(beta_guess)
        
        # Calculate derivative numerically (Central Difference is better, but Forward is fine)
        f_val_plus = theta_beta_m_func(beta_guess + epsilon)
        f_prime = (f_val_plus - f_val) / epsilon
        
        # Newton Step
        beta_new = beta_guess - (f_val / f_prime)
        
        # Error check
        error = abs(beta_new - beta_guess)
        
        # Print step for verification
        print(f"{i:<5} | {math.degrees(beta_new):<12.6f} | {error:<12.6e}")
        
        beta_guess = beta_new
        
        if error < tolerance:
            print("-" * 35)
            print("Converged.")
            return math.degrees(beta_guess)
            
    raise ValueError("Did not converge. Check inputs (Shock might be detached).")

# --- EXECUTION ---
# Case from our discussion
M1 = 36.0
theta_wedge = 15.0
gamma_gas = 1.4

beta_solution = calculate_oblique_shock_beta(M1, theta_wedge, gamma_gas)

print(f"\nFinal Result:")
print(f"Mach: {M1}")
print(f"Wedge Angle: {theta_wedge} deg")
print(f"Shock Angle (Beta): {beta_solution:.4f} deg")

# Verification Calculation of Mn1
beta_rad = math.radians(beta_solution)
Mn1 = M1 * math.sin(beta_rad)
print(f"Normal Mach Component (Mn1): {Mn1:.4f}")