import numpy as np

def gen(
        n: int, 
        g_prob: float = 0.25, 
        fixed: bool = True, 
        dimers: bool = False
    ):
    """
    Randomly generate polymer of length N from L and G monomers
    
    Args:
        N (int): length of polymer (number of monomers it contains)
        g_prob (float): probability of G appearance (default=0.25)
        fixed (bool): whether to fix the number of Gs in string
            - True -> randomly generate g_prob * N fixed positions (no replacement) for G monomers
            - False -> every index has a g_prob probability of being a G 
        dimers (bool): whether polymer should be built by pairs of two polymer dimers 
    
    Returns:
        str: generated polymer
    
    Sample Runs:
        (48, 0.25, True, False) -> LLGLLLGLLLLLGLGLLLLLLLLLLGLLLLLGLGGGGLLGLLLLGLLL
        (48, 0.25, True, True) -> LLLLGGLLLLLLLLLLGGLLGGGGLLLLLLLLLLGGLLLLLLLLLLGG
        (48, 0.25, False, False) -> LLLGGLGLLGLLGLLLLGLLLLLLLLLLLLLGLLGLLLGLLGGGGLLL
    """
    
    if dimers:
        n = int(n / 2)

    polymer_pos = []
        
    if fixed: 
        g_idx = np.random.choice(n, int(n * g_prob), replace=False)

        polymer_pos = np.zeros(n, dtype=int)
        polymer_pos[g_idx] = 1
    else:
        polymer_pos = np.array([1 if np.random.random() <  g_prob else 0 for _ in range(n)])

    polymer = ""
    for i in range(n):
        if not dimers:
            if polymer_pos[i] == 1: polymer += "G"
            else: polymer += "L"
        else:
            if polymer_pos[i] == 1: polymer += "GG"
            else: polymer += "LL"
    
    return polymer

def calc_stats(polymer: str):
    """
    Count number of overlapping GG, LL, GL, LG sub dimers

    Args:
        polymer (str): input polymer of L and G monomers

    Returns:
        4 element tuple: (GGs, LLs, GLs, LGs)
    """
    n = len(polymer)

    num_GGs, num_LLs, num_GLs, num_LGs = 0, 0, 0, 0
    for i in range(n - 1):
        if polymer[i] == 'G' and polymer[i+1] == 'G':
            num_GGs += 1
        if polymer[i] == 'L' and polymer[i+1] == 'L':
            num_LLs += 1
        if polymer[i] == 'G' and polymer[i+1] == 'L':
            num_GLs += 1
        if polymer[i] == 'L' and polymer[i+1] == 'G':
            num_LGs += 1

    return num_GGs, num_LLs, num_GLs, num_LGs

def calc_ratio(polymer: str):
    """
    Calculate the ratio of Gs to Ls in polymer

    Args:
        polymer (str): input polymer of L and G monomers

    Returns:
        float: ratio of Gs to Ls
    """
    num_Gs = polymer.count('G')
    num_Ls = polymer.count('L')
    
    return num_Gs / num_Ls

def main():
    N = 10000
    p_length = 100
    
    GGs, LLs, GLs, LGs = [], [], [], []

    ratios = []
    with open(f'data/sample_polymers_{p_length}.out', 'w') as f:
        for _ in range(N):
            polymer = gen(p_length, fixed=True, dimers=False)
            f.write(polymer + '\n')

            (GG, LL, GL, LG) = calc_stats(polymer)
            GGs.append(GG)
            LLs.append(LL)
            GLs.append(GL)
            LGs.append(float(LG))

            ratios.append(calc_ratio(polymer))

    GGs = np.array(GGs)
    LLs = np.array(LLs)
    GLs = np.array(GLs)
    LGs = np.array(LGs)

    print(f"n: {p_length}")
    print("-" * 20)
    print(f"Mean G-Gs: {np.mean(GGs):.2f}")
    print(f"SEM G-Gs:  {np.std(GGs) / np.sqrt(N - 1):.2f}")
    print("-" * 20)

    print(f"Mean L-Ls: {np.mean(LLs):.2f}")
    print(f"SEM L-Ls:  {np.std(LLs) / np.sqrt(N - 1):.2f}")
    print("-" * 20)

    print(f"Mean G-Ls: {np.mean(GLs):.2f}")
    print(f"SEM G-Ls:  {np.std(GLs) / np.sqrt(N - 1):.2f}")
    print("-" * 20)

    print(f"Mean L-Gs: {np.mean(LGs):.2f}")
    print(f"SEM L-Gs:  {np.std(LGs) / np.sqrt(N - 1):.2f}")
    print("-" * 20)
    
    # Correct 0s to 1s to avoid divide by 0
    for i in range(len(LGs)):
        if(LGs[i] == 0):
            LGs[i] = 1
        if(GLs[i] == 0):
            GLs[i] = 1
    
    print("OTHER METRICS")
    print(f"L_L (mean)    = {np.mean(LLs / LGs + 1):.2f}")
    print(f"L_L (sem)     = {np.std(LLs / LGs + 1) / np.sqrt(N - 1):.2f}")
    print()
    print(f"L_G (mean)    = {np.mean(GGs / GLs + 1):.2f}")
    print(f"L_G (sem)     = {np.std(GGs / GLs + 1) / np.sqrt(N - 1):.2f}")
    print()
    print(f"R_c (mean)    = {np.mean(GGs / GLs):.2f}")
    print(f"R_c (sem)     = {np.std(GGs / GLs) / np.sqrt(N - 1):.2f}")
    print()
    print(f"Ratio of L/G (mean) = {np.mean(ratios):.2f}")

if __name__ == "__main__":
    main()