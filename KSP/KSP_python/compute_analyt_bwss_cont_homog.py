import numpy as np
import matplotlib.pyplot as plt
import random
import math
import argparse


# assumes linear dependence on gamma
def integrand(r_s, th_s, r_t, th_t, r_v, th_v, gamma, eps):
    
    sx = r_s * np.cos(th_s)
    sy = r_s * np.sin(th_s)
    tx = r_t * np.cos(th_t)
    ty = r_t * np.sin(th_t)
    vx = r_v * np.cos(th_v)
    vy = r_v * np.sin(th_v)

    if eps > 0:
        d_st = np.sqrt((sx - tx)**2 + (sy - ty)**2)
        d_sv = np.sqrt((vx - sx)**2 + (vy - sy)**2)
        d_tv = np.sqrt((vx - tx)**2 + (vy - ty)**2)

        diff = d_st * (1 + eps) - d_sv - d_tv
        norm = np.exp(gamma * eps * d_st)

        # Apply the conditional formula: only for diff ≥ 0
        positive_mask = diff >= 0
        result = np.zeros_like(diff)
        gd = gamma * diff[positive_mask]
        result[positive_mask] = (np.exp(gd) * (gd + 1)) / norm[positive_mask]
        return result
    
    else: # Shortest path case (eps = 0)
        B = 2 * (sy * (vy - sy) + sx * (vx - sx))
        A = (vx - sx)**2 + (vy - sy)**2
        C = sx**2 + sy**2 - 1
        discriminant = B**2 - 4*A*C

        result = np.zeros_like(discriminant)

        real_mask = discriminant >= 0
        sqrt_disc = np.zeros_like(discriminant)
        sqrt_disc[real_mask] = np.sqrt(discriminant[real_mask])

        lambda_max = np.zeros_like(discriminant)
        lambda_max[real_mask] = (-B[real_mask] + sqrt_disc[real_mask]) / (2 * A[real_mask])

        # Integrand is max(lambda_max - 1, 0)
        result[real_mask] = np.maximum(lambda_max[real_mask] - 1, 0)
        return result



def monte_carlo_integrate_polar(f, samples=100000, args=()):
    # Sample u uniformly in [0,1], then r = sqrt(u) to match area element
    u_rs = np.random.rand(samples)
    u_rt = np.random.rand(samples)
    r_s = np.sqrt(u_rs)
    r_t = np.sqrt(u_rt)
    
    th_s = np.random.uniform(0, 2*np.pi, samples)
    th_t = np.random.uniform(0, 2*np.pi, samples)
    eps = args[-1] if len(args) > 0 else 0
    # Compute the integrand values
    # values = np.array([
    #     f(r_s[i], th_s[i], r_t[i], th_t[i], *args) 
    #     for i in range(samples)
    # ])
    values = f(r_s, th_s, r_t, th_t, *args)

    
    # The volume of integration is:
    if eps > 0:
        volume = (np.pi)**2  # for r_s^2 /2 and theta_s
    elif eps == 0:
         volume = np.pi
    estimate = volume * np.mean(values)
    std_error = volume * np.std(values) / np.sqrt(samples)
    
    return estimate, std_error


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Compute analytical betweenness centrality for homogeneous continuous networks.")
    parser.add_argument("gamma", type=float, help="Related to the increase of paths with distance")
    parser.add_argument("eps", type=float, help="Tolerance of the QSP-BW")

    args = parser.parse_args()

    r = np.linspace(0, 1, 60)         
    theta = np.linspace(0, 2*np.pi, 50) 
    # Create a meshgrid for r and theta
    bwss = np.zeros((len(r), len(theta)))
    errs = np.zeros((len(r), len(theta)))

    for i, r_v in enumerate(r):
        print(f"Integrating for r_v = {r_v:.2f}...")
        for j, th_v in enumerate(theta):

            result, err = monte_carlo_integrate_polar(integrand, args=[r_v, th_v, args.gamma, args.eps])
            bwss[i, j] = result
            errs[i, j] = err


    # save to file
    # np.savetxt(f'./bwss_continuum/bwss_polar_{int(args.eps*100)}_gam{int(args.gamma)}.txt', bwss, delimiter=',')
    # np.savetxt(f'./bwss_continuum/err_bwss_polar_{int(args.eps*100)}_gam{int(args.gamma)}.txt', errs, delimiter=',')

    np.savetxt(f'./bwss_continuum/bwss_polar_{int(args.eps*100)}_gam{int(args.gamma)}.txt', bwss, delimiter=',')
    np.savetxt(f'./bwss_continuum/err_bwss_polar_{int(args.eps*100)}_gam{int(args.gamma)}.txt', errs, delimiter=',')

    print("Integration complete. Results saved to files.")

    # def integrand(r_s, th_s, r_t, th_t, r_v, th_v, gamma, eps):
#     # Convert polar to cartesian
#     sx = r_s * np.cos(th_s)
#     sy = r_s * np.sin(th_s)
#     tx = r_t * np.cos(th_t)
#     ty = r_t * np.sin(th_t)
#     vx = r_v * np.cos(th_v)
#     vy = r_v * np.sin(th_v)

#     if eps > 0:
#         d_st = np.sqrt((sx - tx)**2 + (sy - ty)**2)
#         d_sv = np.sqrt((vx - sx)**2 + (vy - sy)**2)
#         d_tv = np.sqrt((vx - tx)**2 + (vy - ty)**2)

#         diff = d_st * (1 + eps) - d_sv - d_tv
#         norm = np.exp(gamma * eps * d_st)

#         # Apply the conditional formula: only for diff ≥ 0
#         positive_mask = diff >= 0
#         result = np.zeros_like(diff)
#         gd = gamma * diff[positive_mask]
#         result[positive_mask] = (np.exp(gd) * (gd - 1) + 1) / norm[positive_mask]
#         return result

#     else: # Shortest path case (eps = 0)
#         B = 2 * (sy * (vy - sy) + sx * (vx - sx))
#         A = (vx - sx)**2 + (vy - sy)**2
#         C = sx**2 + sy**2 - 1
#         discriminant = B**2 - 4*A*C

#         result = np.zeros_like(discriminant)

#         real_mask = discriminant >= 0
#         sqrt_disc = np.zeros_like(discriminant)
#         sqrt_disc[real_mask] = np.sqrt(discriminant[real_mask])

#         lambda_max = np.zeros_like(discriminant)
#         lambda_max[real_mask] = (-B[real_mask] + sqrt_disc[real_mask]) / (2 * A[real_mask])

#         # Integrand is max(lambda_max - 1, 0)
#         result[real_mask] = np.maximum(lambda_max[real_mask] - 1, 0)
#         return result
    