#!/usr/bin/env python3
"""
Output statistics from vcf
VCF class and some summary statistic functions originate from BaSe (https://github.com/ulasisik/balancing-selection)
Modifications made for use with SWIF(r) by K. Ahlquist
Modifications from N. Font-Porterias, Feb 2025: 
- To increase speed: np.arange() instead of range(), _crop_for_target(), mean_pwisedis, fusion function calc_r2_and_kelly_zns() LD+Zns, 
- Other modifications: calc_zeng_e calculation fix, add a new argument for window size, default value for n argument for sample size in VCF() to None (instead of 17),  
"""
import os
import sys
import allel
import random
import numpy as np
import pandas as pd
import time
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors


def calc_median_r2(g):
    """Calculates median LD r^2"""
    gn = g.to_n_alt(fill=-1)
    LDr = allel.rogers_huff_r(gn)
    LDr2 = LDr ** 2
    median_r2 = np.nanmedian(LDr2)
    return median_r2

def calc_kelly_zns(g, n_pos):
    """Calculates Kelly's Zns statistic"""
    gn = g.to_n_alt(fill=-1)
    LDr = allel.rogers_huff_r(gn)
    LDr2 = LDr ** 2
    kellyzn = (np.nansum(LDr2) * 2.0) / (n_pos * (n_pos - 1.0))
    return kellyzn

def calc_r2_and_kelly_zns(g, n_pos):
    """Calculates median LD r^2 and Kelly's Zns statistic"""
    gn = g.to_n_alt(fill=-1)
    LDr = allel.rogers_huff_r(gn)
    LDr2 = LDr ** 2

    # Calculate median r^2
    median_r2 = np.nanmedian(LDr2)

    # Calculate Kelly's Zns
    kellyzn = (np.nansum(LDr2) * 2.0) / (n_pos * (n_pos - 1.0))

    return median_r2, kellyzn

def calc_pi(croms):
    """Calculates pi"""
    dis1 = []
    for i in np.arange(croms.shape[0]):
        d1 = []
        for j in np.arange(i + 1, croms.shape[0]):
            d1.append(sum(croms[i, :] != croms[j, :]))
        dis1.append(sum(d1))
    pi_est1 = (sum(dis1) / ((croms.shape[0] * (croms.shape[0] - 1.0)) / 2.0))
    return pi_est1

def calc_faywu_h(croms):
    """Calculates Fay and Wu's H statistic"""
    n_sam1 = croms.shape[0]
    counts1 = croms.sum(axis=0)
    S_i1 = []
    for i in np.arange(1, n_sam1):
        S_i1.append(sum(counts1 == i))
    i = np.arange(1, n_sam1)
    n_i = np.subtract(n_sam1, i)
    thetaP1 = sum((n_i * i * S_i1 * 2) / (n_sam1 * (n_sam1 - 1.0)))
    thetaH1 = sum((2 * np.multiply(S_i1, np.power(i, 2))) / (n_sam1 * (n_sam1 - 1.0)))
    Hstat1 = thetaP1 - thetaH1
    return Hstat1

def calc_fuli_f_star(croms):
    """Calculates Fu and Li's D* statistic"""
    n_sam1 = croms.shape[0]
    n_pos1 = np.size(croms, 1)
    an = np.sum(np.divide(1.0, np.arange(1, n_sam1)))
    bn = np.sum(np.divide(1.0, np.power(np.arange(1, n_sam1), 2)))
    an1 = an + np.true_divide(1, n_sam1)

    vfs = (((2 * (n_sam1 ** 3.0) + 110.0 * (n_sam1 ** 2.0) - 255.0 * n_sam1 + 153) / (
            9 * (n_sam1 ** 2.0) * (n_sam1 - 1.0))) + ((2 * (n_sam1 - 1.0) * an) / (n_sam1 ** 2.0)) - (
                   (8.0 * bn) / n_sam1)) / ((an ** 2.0) + bn)
    ufs = ((n_sam1 / (n_sam1 + 1.0) + (n_sam1 + 1.0) / (3 * (n_sam1 - 1.0)) - 4.0 / (
            n_sam1 * (n_sam1 - 1.0)) + ((2 * (n_sam1 + 1.0)) / ((n_sam1 - 1.0) ** 2)) * (
                    an1 - ((2.0 * n_sam1) / (n_sam1 + 1.0)))) / an) - vfs

    pi_est = calc_pi(croms)
    ss = sum(np.sum(croms, axis=0) == 1)
    Fstar1 = (pi_est - (((n_sam1 - 1.0) / n_sam1) * ss)) / ((ufs * n_pos1 + vfs * (n_pos1 ** 2.0)) ** 0.5)
    return Fstar1

def calc_fuli_d_star(croms):
    """Calculates Fu and Li's D* statistic"""
    n_sam1 = croms.shape[0]
    n_pos1 = np.size(croms, 1)
    an = np.sum(np.divide(1.0, np.arange(1, n_sam1)))
    bn = np.sum(np.divide(1.0, np.power(np.arange(1, n_sam1), 2)))
    an1 = an + np.true_divide(1, n_sam1)

    cn = (2 * (((n_sam1 * an) - 2 * (n_sam1 - 1))) / ((n_sam1 - 1) * (n_sam1 - 2)))
    dn = (cn + np.true_divide((n_sam1 - 2), ((n_sam1 - 1) ** 2)) + np.true_divide(2, (n_sam1 - 1)) * (
            3.0 / 2 - (2 * an1 - 3) / (n_sam1 - 2) - 1.0 / n_sam1))

    vds = (((n_sam1 / (n_sam1 - 1.0)) ** 2) * bn + (an ** 2) * dn - 2 * (n_sam1 * an * (an + 1)) / (
            (n_sam1 - 1.0) ** 2)) / (an ** 2 + bn)
    uds = ((n_sam1 / (n_sam1 - 1.0)) * (an - n_sam1 / (n_sam1 - 1.0))) - vds

    ss = sum(np.sum(croms, axis=0) == 1)
    Dstar1 = ((n_sam1 / (n_sam1 - 1.0)) * n_pos1 - (an * ss)) / (uds * n_pos1 + vds * (n_pos1 ^ 2)) ** 0.5
    return Dstar1

def calc_zeng_e(croms):
    """Calculates Zeng et al's E statistic"""
    n_sam1 = croms.shape[0]
    n_pos1 = np.size(croms, 1)
    an = np.sum(np.divide(1.0, np.arange(1, n_sam1)))
    bn = np.sum(np.divide(1.0, np.power(np.arange(1, n_sam1), 2)))
    counts1 = croms.sum(axis=0)

    S = sum((counts1 > 0) & (counts1 < n_sam1))# Count the number of segregating sites (where at least one sample has a different genotype)
    thetaW = S / an

    S_i1 = []
    for i in np.arange(1, n_sam1):
        S_i1.append(sum(counts1 == i))
    thetaL = np.sum(np.multiply(S_i1, np.arange(1, n_sam1))) / (n_sam1 - 1.0)
    theta2 = (S * (S - 1.0)) / (an ** 2 + bn)

    var1 = (n_sam1 / (2.0 * (n_sam1 - 1.0)) - 1.0 / an) * thetaW
    var2 = theta2 * (bn / (an ** 2.0)) + 2 * bn * (n_sam1 / (n_sam1 - 1.0)) ** 2.0 - (2.0 * (n_sam1 * bn - n_sam1 + 1.0)) / ((n_sam1 - 1.0) * an) - (3.0 * n_sam1 + 1.0) / (n_sam1 - 1.0)
    varlw = var1 + var2

    ZengE1 = (thetaL - thetaW) / (varlw) ** 0.5
    return ZengE1

def calc_rageddness(croms):
    """Calculates rageddness index"""
    mist = []
    for i in np.arange(croms.shape[0] - 1):
        for j in np.arange(i + 1, croms.shape[0]):
            mist.append(sum(croms[i, :] != croms[j, :]))
    mist = np.array(mist)
    lnt = mist.shape[0]
    fclass = []
    for i in np.arange(1, np.max(mist) + 2):
        fclass.append((np.true_divide(sum(mist == i), lnt) - np.true_divide(sum(mist == (i - 1)), lnt)) ** 2)
    rgd1 = np.sum(fclass)
    return rgd1

def sum_stats(croms, pos, NCHROMS, target_pos, file_name, croms_exact):
    """
    Calculates summary statistics
    Parameters:
        croms: haplotype matrix
        pos: positions
        sname: simulation name
        NCHROMS: number of chromosomes
    Returns:
        A list of labels (names of statistics)
        A list of values
    """

    # SUMMARY STATISTICS
    # REGION 1: central 10kb([20kb:30kb])
    #pos1 = pos[np.logical_and(pos > 20000, pos < 30000)]
    #croms1 = croms[:, np.logical_and(pos > 20000, pos < 30000).tolist()]

    # Window size = 100kb (instead of 10kb)
    pos1 = pos[np.logical_and(pos > 10000, pos < 110000)]
    croms1 = croms[:, np.logical_and(pos > 10000, pos < 110000).tolist()]

    n_pos1 = np.size(croms1, 1)
    freq1 = np.true_divide(croms1.sum(axis=0), NCHROMS)
    freq1 = np.array(freq1)

    haplos = np.transpose(croms1)
    haplos_exact = np.transpose(croms_exact)

    h1 = allel.HaplotypeArray(haplos)
    h_e = allel.HaplotypeArray(haplos_exact)

    ac1 = h1.count_alleles()
    ac_e = h_e.count_alleles()

    g1 = h1.to_genotypes(ploidy=2, copy=True)

    g_e = h_e.to_genotypes(ploidy=2, copy=True)


    # mean_pairwise_distance
    pwise_dis = allel.mean_pairwise_difference(ac_e)[0]
    mean_pwise_dis1 = allel.mean_pairwise_difference(ac1)
    mean_mean_pwise_dis1, median_mean_pwise_dis1, max_mean_pwise_dis1 = (
        np.mean(mean_pwise_dis1),
        np.median(mean_pwise_dis1),
        np.max(mean_pwise_dis1)
    )

    # tajimasd
    TjD1 = allel.tajima_d(ac1)

    # watterson
    theta_hat_w1 = allel.watterson_theta(pos1, ac1)

    # heterogeneity
    obs_het_site = allel.heterozygosity_observed(g_e)[0]
    obs_het1 = allel.heterozygosity_observed(g1)
    af1 = ac1.to_frequencies()
    exp_het1 = allel.heterozygosity_expected(af1, ploidy=2)
    mean_obs_het1 = np.mean(obs_het1)
    median_obs_het1 = np.median(obs_het1)
    max_obs_het1 = np.max(obs_het1)
    # return 0, if 0/0 encountered
    ob_exp_het1 = np.divide(obs_het1, exp_het1, out=np.zeros_like(obs_het1), where=exp_het1 != 0)
    mean_obs_exp1 = np.nanmean(ob_exp_het1)
    median_obs_exp1 = np.nanmedian(ob_exp_het1)
    max_obs_exp1 = np.nanmax(ob_exp_het1)

    # Haplotype_stats
    hh1 = allel.garud_h(h1)
    h11 = hh1[0]
    h121 = hh1[1]
    h1231 = hh1[2]
    h2_h11 = hh1[3]
    n_hap1 = np.unique(croms1, axis=0).shape[0]
    hap_div1 = allel.haplotype_diversity(h1)

    ehh1 = allel.ehh_decay(h1)
    mean_ehh1 = np.mean(ehh1)
    median_ehh1 = np.median(ehh1)

    ihs1 = allel.ihs(h1, pos1, include_edges=True)
    median_ihs1 = np.nanmedian(ihs1)

    # nsl
    nsl1 = allel.nsl(h1)
    max_nsl1 = np.nanmax(nsl1)
    median_nsl1 = np.nanmedian(nsl1)

    # NCD
    tf = 0.5
    freq11 = freq1[freq1 < 1]
    n1 = freq11.shape[0]
    ncd11 = (sum((freq11 - tf) ** 2) / n1) ** 0.5

    # kelly Zns and median LD r2
    median_r21, kellyzn1 = calc_r2_and_kelly_zns(g1, n_pos1)

    # pi
    #pi_est1 = calc_pi(croms1)

    # FayWusH
    Hstat1 = calc_faywu_h(croms1)

    # of singletons
    Ss1 = sum(np.sum(croms1, axis=0) == 1)

    # fu_li Dstar
    #Dstar1 = calc_fuli_d_star(croms1)

    # fu_li Fstar
    #Fstar1 = calc_fuli_f_star(croms1)

    # Zeng_E
    ZengE1 = calc_zeng_e(croms1)

    # rageddness index
    #rgd1 = calc_rageddness(croms1)

    stats = [target_pos, file_name, pwise_dis, mean_mean_pwise_dis1, median_mean_pwise_dis1, max_mean_pwise_dis1,
             TjD1, theta_hat_w1, obs_het_site, mean_obs_het1, median_obs_het1, max_obs_het1, mean_obs_exp1, median_obs_exp1,
             max_obs_exp1, median_r21, h11, h121, h1231, h2_h11, hap_div1, n_hap1, mean_ehh1, median_ehh1, median_ihs1,
             max_nsl1, median_nsl1, ncd11, kellyzn1, Hstat1, Ss1, ZengE1]

    labs = ['Position', 'FileName', 'PwiseDist', 'Mean(MeanPwiseDist)1', 'Median(MeanPwiseDist)1', 'Max(MeanPwiseDist)1',
            'Tajimas D1', 'Watterson1', 'ObservedHet', 'Mean(ObservedHet)1', 'Median(ObservedHet)1', 'Max(ObservedHet)1',
            'Mean(Obs/Exp Het)1', 'Median(Obs/Exp Het)1', 'Max(Obs/Exp Het)1', 'Median(r2)1', 'H1_1', 'H12_1',
            'H123_1', 'H2/H1_1', 'Haplotype Diversity1', '# of Hap1', 'Mean(EHH)1', 'Median(EHH)1', 'Median(ihs)1',
            'Max(nsl)1', 'Median(nsl)1', 'NCD1_1', 'KellyZns1', 'faywuH1', '#ofSingletons1', 'ZengE1']

    return labs, stats


class VCF(object):

    def __init__(self, file_name, n=None, seed=None):
        """Loads VCF file
        Parameters:
            file_name: target vcf file
            n: number of diploid individuals to be selected randomly
            seed: random seed
        """
        startTime = time.time()

        self.file_name=file_name
        vcfile = allel.read_vcf(file_name)

        positions = vcfile["variants/POS"].tolist()
        snp_id = vcfile["variants/POS"].tolist() #Since ID is nonspecific in my data, changed ID to POS

        gts = allel.GenotypeArray(vcfile["calldata/GT"])

        # Randomly select n samples
        if n != gts.shape[1]:
            random.seed(seed)
            samp_to_select = random.sample(range(0, gts.shape[1]), k=int(n))

            gts = gts.subset(sel0=None, sel1=samp_to_select)

        hps = gts.to_haplotypes()
        haplotypes = np.transpose(hps.values)

        # Filter non biallelic positions
        nonbiallel_pos = np.sum(haplotypes >= 2, axis=0) == 0
        if len(nonbiallel_pos) > sum(nonbiallel_pos):
            print("Filtering {} non-biallelic positions..".format(len(nonbiallel_pos) - sum(nonbiallel_pos)))
            haplotypes = haplotypes[:, nonbiallel_pos]
            positions = np.array(positions)[nonbiallel_pos].tolist()
            snp_id = np.array(snp_id)[nonbiallel_pos].tolist()

        if gts.count_missing():
            print("Warning: {} missing positions are set to 0 (ancestral)".format(gts.count_missing()))
            mismask = haplotypes < 0
            haplotypes[mismask] = 0
        
        print(f'Time to load VCF: {time.time() - startTime}')

        self.croms, self.positions, self.snp_id = haplotypes, positions, snp_id

    def scan_targets(self, N, target_range, target_freq):
        """Scans vcf file for candidate targets
        Parameters:
            N: length of the (simulated) sequence
            target_range: A tuple specifying the target range of positions. If None, scans all the positions.
            target_freq: A tuple specifying the frequency range for targets.
        Returns:
            A list of candidate targets
        """
        if not target_range:
            target_range = [int(min(self.positions) + N // 2), int(max(self.positions) - N // 2)]
        else:
            target_range = [int(target_range[0] + N // 2), int(target_range[1] - N // 2)] #add also "buffer" to avoid error "Flanking regions..."
        print("target_range_new")    
        print(target_range)
        freqs = np.true_divide(np.sum(self.croms, axis=0), self.croms.shape[0])

        targets = np.logical_and(np.logical_and(freqs > target_freq[0], freqs < target_freq[1]),
                                 np.logical_and(np.array(self.positions) >= target_range[0],
                                                np.array(self.positions) <= target_range[1]))
        if len(targets) == 0:
            raise LookupError("No SNP found for the given target frequency range")

        target_list = np.array(self.snp_id)[targets]

        return target_list

    def _crop_for_target(self, N, target_snp):
        """
        Given the target SNP, crops the sequence so that the target SNP is at the center and length of the sequence
        is equal to N.
        Parameters:
            N: length of the (simulated) sequence
            target_snp: target SNP id
        Return:
            Haplotype matrix
            Positions of segregating sites
            SNP IDs
        """

        # Convert positions to numpy array once
        positions = np.array(self.positions)
        snp_idx = self.snp_id.index(target_snp)  # index
        pos_mut = positions[snp_idx]  # actual VCF position from index

        up_bound = pos_mut + N / 2
        low_bound = pos_mut - N / 2

        if low_bound < 0 or up_bound > positions[-1]:  # Assuming positions is sorted
            raise IndexError(f"Flanking regions around the target SNP must be greater than {N // 2}")

        # Use np.searchsorted to find indices directly
        low_idx = np.searchsorted(positions, low_bound)
        up_idx = np.searchsorted(positions, up_bound)

        # Get the new positions and haplotypes
        pos_new = positions[low_idx:up_idx] + N / 2 - pos_mut
        croms_new = self.croms[:, low_idx:up_idx]

        # Get GT values for the target SNP
        exact_idx = np.where(positions == pos_mut)[0]
        if exact_idx.size > 0:
            croms_exact = self.croms[:, exact_idx][:, :1]
        else:
            croms_exact = np.empty((self.croms.shape[0], 0))

        new_snp_id = np.array(self.snp_id)[low_idx:up_idx]

        return croms_new, pos_new, new_snp_id, croms_exact

    def create_stat(self, N, target_freq=(0.4, 0.6), target_list=None, target_range=None, scale=False, pca=False, output=True, output_name=None):
        """
        Scans the VCF file and create and preprocess summary statistics
        Parameters:
            N: length of the (simulated) sequence
            target_freq: A tuple specifying the frequency range for targets
            target_list: A list of target SNPs. If None, scans the target region for all candidate targets
            target_range:  A tuple specifying the target range of positions. If None, scans all the positions
            scale: If True, performs feature scaling
            pca: If True, performs pca
            w: window size (default 50kb)

        Returns:
            A matrix containing summary statistics
            A list containing corresponding SNP IDs
        """
        if not target_list:
            startTime = time.time()
            target_list = self.scan_targets(N=N, target_range=target_range, target_freq=target_freq)
            print("{} candidate targets have been found.".format(len(target_list)))
            print(f'Time to scan target: {time.time() - startTime}')

        else:
            if not all(t in self.snp_id for t in target_list):
                raise ValueError("Target SNP not found in VCF file")

        file = self.file_name.split("/")[-1]
        statdf = pd.DataFrame()
        target_pos = []
        for target_snp in target_list:
            startTime = time.time()
            print("Processing position: " + str(target_snp))
            croms, pos, _ , croms_exact= self._crop_for_target(N, target_snp)
            labs, stats = sum_stats(croms, pos, croms.shape[0], str(target_snp), file, croms_exact)
            statdf = statdf.append(dict(zip(labs, stats)), ignore_index=True)
            statdf = statdf[labs]
            target_pos.append(self.positions[self.snp_id.index(target_snp)])
            posTime = (time.time() - startTime)
            print('Time to process position: ' + str(posTime))
        stat_matrix = statdf.iloc[:, 3:].values

        if scale:
            sc = StandardScaler()
            sc.fit(stat_matrix)
            stat_matrix = sc.transform(stat_matrix)

        if pca:
            pca = PCA()
            pca.fit(stat_matrix)
            stat_matrix_pca = pca.transform(stat_matrix)
            stat_matrix = stat_matrix_pca
        
        if output:
            if output_name:
                print("Outputting to file: " + output_name)
                statdf.to_csv(output_name)
            else:
                print("Outputting to file: " + self.file_name[:-3] + '_unnormed_stats.csv')
                statdf.to_csv(self.file_name[:-3] + '_unnormed_stats.csv')

        return stat_matrix, target_list, target_pos
