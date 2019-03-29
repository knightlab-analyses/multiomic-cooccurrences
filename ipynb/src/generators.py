import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.utils import check_random_state
from skbio.stats.composition import clr_inv as softmax
from skbio.stats.composition import ilr_inv, clr, closure
from biom.util import biom_open
from biom import Table
from .sim import partition_microbes, partition_metabolites
from scipy.stats import norm, expon


# utility functions for saving files
def deposit_biofilms(output_dir,
                     abs_table1, abs_table2,
                     rel_table1, rel_table2,
                     edges, metadata, sample_id):
    """ Writes down tables and edges into files.

    Parameters
    ----------
    output_dir : str
        output directory
    rel_table1 : biom.Table
        Biom table of relative abundances
    rel_table2 : biom.Table
        Biom table of relative abundances
    abs_table1 : biom.Table
        Biom table of absolute abundances
    abs_table2 : biom.Table
        Biom table of absolute abundances
    edges : list
        Edge list for ground truthing.
    metadata : pd.DataFrame
        Dataframe of sample metadata
    sample_id : str
        sample id
    """
    output_abs_microbes = "%s/table.abs.microbes.%s.biom" % (
        output_dir, sample_id)
    output_abs_metabolites = "%s/table.abs.metabolites.%s.biom" % (
        output_dir, sample_id)
    output_rel_microbes = "%s/table.rel.microbes.%s.biom" % (
        output_dir, sample_id)
    output_rel_metabolites = "%s/table.rel.metabolites.%s.biom" % (
        output_dir, sample_id)
    output_md = "%s/metadata.%s.txt" % (
        output_dir, sample_id)
    output_U = "%s/U.%s.txt" % (
        output_dir, sample_id)
    output_V = "%s/V.%s.txt" % (
        output_dir, sample_id)
    output_edges = "%s/edges.%s.txt" % (
        output_dir, sample_id)
    output_ranks = "%s/ranks.%s.txt" % (
        output_dir, sample_id)

    # idx1 = table1.sum(axis=0) > 0
    # idx2 = table2.sum(axis=0) > 0
    # table1 = table1.loc[:, idx1]
    # table2 = table2.loc[:, idx2]

    # relative abundances
    table1 = Table(rel_table1.values.T, rel_table1.columns, rel_table1.index)
    table2 = Table(rel_table2.values.T, rel_table2.columns, rel_table2.index)
    with biom_open(output_rel_microbes, 'w') as f:
        table1.to_hdf5(f, generated_by='moi1')
    with biom_open(output_rel_metabolites, 'w') as f:
        table2.to_hdf5(f, generated_by='moi2')

    # absolute abundances
    table1 = Table(abs_table1.values.T, abs_table1.columns, abs_table1.index)
    table2 = Table(abs_table2.values.T, abs_table2.columns, abs_table2.index)
    with biom_open(output_abs_microbes, 'w') as f:
        table1.to_hdf5(f, generated_by='moi1')
    with biom_open(output_abs_metabolites, 'w') as f:
        table2.to_hdf5(f, generated_by='moi2')

    pd.DataFrame(edges).to_csv(output_edges, sep='\t')
    metadata.to_csv(output_md, sep='\t')


def deposit(output_dir, table1, table2, metadata, U, V, B, it, rep):
    """ Writes down tables, metadata and feature metadata into files.

    Parameters
    ----------
    output_dir : str
        output directory
    table1 : biom.Table
        Biom table
    table2 : biom.Table
        Biom table
    metadata : pd.DataFrame
        Dataframe of sample metadata
    U : np.array
        Microbial latent variables
    V : np.array
        Metabolite latent variables
    edges : list
        Edge list for ground truthing.
    feature_metadata : pd.DataFrame
        Dataframe of features metadata
    it : int
        iteration number
    rep : int
        repetition number
    """
    choice = 'abcdefghijklmnopqrstuvwxyz'
    output_microbes = "%s/table_microbes.%d_%s.biom" % (
        output_dir, it, choice[rep])
    output_metabolites = "%s/table_metabolites.%d_%s.biom" % (
        output_dir, it, choice[rep])
    output_md = "%s/metadata.%d_%s.txt" % (
        output_dir, it, choice[rep])
    output_U = "%s/U.%d_%s.txt" % (
        output_dir, it, choice[rep])
    output_V = "%s/V.%d_%s.txt" % (
        output_dir, it, choice[rep])
    output_B = "%s/B.%d_%s.txt" % (
        output_dir, it, choice[rep])
    output_ranks = "%s/ranks.%d_%s.txt" % (
        output_dir, it, choice[rep])

    idx1 = table1.sum(axis=0) > 0
    idx2 = table2.sum(axis=0) > 0
    table1 = table1.loc[:, idx1]
    table2 = table2.loc[:, idx2]

    table1 = Table(table1.values.T, table1.columns, table1.index)
    table2 = Table(table2.values.T, table2.columns, table2.index)

    with biom_open(output_microbes, 'w') as f:
        table1.to_hdf5(f, generated_by='moi1')
    with biom_open(output_metabolites, 'w') as f:
        table2.to_hdf5(f, generated_by='moi2')

    ranks = clr(softmax(np.hstack(
        (np.zeros((U.shape[0], 1)), U @ V))))
    ranks = ranks[idx1, :]
    ranks = ranks[:, idx2]
    ranks = pd.DataFrame(
        ranks, index=table1.ids(axis='observation'),
        columns=table2.ids(axis='observation'))
    ranks.to_csv(output_ranks, sep='\t')
    metadata.to_csv(output_md, sep='\t', index_label='#SampleID')

    np.savetxt(output_B, B)
    np.savetxt(output_U, U)
    np.savetxt(output_V, V)


def deposit_biofilm(table1, table2, metadata, U, V, edges, it, rep, output_dir):
    """ Writes down tables, metadata and feature metadata into files.

    Parameters
    ----------
    table : biom.Table
        Biom table
    metadata : pd.DataFrame
        Dataframe of sample metadata
    feature_metadata : pd.DataFrame
        Dataframe of features metadata
    it : int
        iteration number
    rep : int
        repetition number
    output_dir : str
        output directory
    """
    choice = 'abcdefghijklmnopqrstuvwxyz'
    output_microbes = "%s/table_microbes.%d_%s.biom" % (
        output_dir, it, choice[rep])
    output_metabolites = "%s/table_metabolites.%d_%s.biom" % (
        output_dir, it, choice[rep])
    output_md = "%s/metadata.%d_%s.txt" % (
        output_dir, it, choice[rep])
    output_U = "%s/U.%d_%s.txt" % (
        output_dir, it, choice[rep])
    output_V = "%s/V.%d_%s.txt" % (
        output_dir, it, choice[rep])
    output_B = "%s/edges.%d_%s.txt" % (
        output_dir, it, choice[rep])
    output_ranks = "%s/ranks.%d_%s.txt" % (
        output_dir, it, choice[rep])

    idx1 = table1.sum(axis=0) > 0
    idx2 = table2.sum(axis=0) > 0
    table1 = table1.loc[:, idx1]
    table2 = table2.loc[:, idx2]

    table1 = Table(table1.values.T, table1.columns, table1.index)
    table2 = Table(table2.values.T, table2.columns, table2.index)

    with biom_open(output_microbes, 'w') as f:
        table1.to_hdf5(f, generated_by='moi1')
    with biom_open(output_metabolites, 'w') as f:
        table2.to_hdf5(f, generated_by='moi2')

    ranks = (U @ V)

    ranks = ranks[idx1, :]
    ranks = ranks[:, idx2]
    ranks = pd.DataFrame(ranks, index=table1.ids(axis='observation'),
                         columns=table2.ids(axis='observation'))
    ranks.to_csv(output_ranks, sep='\t')
    metadata.to_csv(output_md, sep='\t', index_label='#SampleID')

    B = B[:, idx1]

    np.savetxt(output_U, U)
    np.savetxt(output_V, V)
    np.savetxt(output_B, B)


def deposit_blocktable(output_dir, abs_table, rel_table, metadata, truth, sample_id):
    choice = 'abcdefghijklmnopqrstuvwxyz'
    output_abstable = "%s/rel_table.%s.biom" % (
        output_dir, sample_id)
    output_reltable = "%s/abs_table.%s.biom" % (
        output_dir, sample_id)
    output_metadata = "%s/metadata.%s.txt" % (
        output_dir, sample_id)
    output_truth = "%s/truth.%s.txt" % (
        output_dir, sample_id)

    abs_t = Table(abs_table.T.values,
                  abs_table.columns.values,
                  abs_table.index.values)
    with biom_open(output_abstable, 'w') as f:
        abs_t.to_hdf5(f, generated_by='moi')

    rel_t = Table(rel_table.T.values,
                  rel_table.columns.values,
                  rel_table.index.values)
    with biom_open(output_reltable, 'w') as f:
        rel_t.to_hdf5(f, generated_by='moi')

    metadata.to_csv(output_metadata, sep='\t')
    truth.to_csv(output_truth, sep='\t')


# This is for multiomics benchmarks
def random_sigmoid_multimodal(
        num_microbes=20, num_metabolites=100, num_samples=100,
        num_latent_microbes=5, num_latent_metabolites=10,
        num_latent_shared=3, low=-1, high=1,
        microbe_total=10, metabolite_total=100,
        uB=0, sigmaB=2, sigmaQ=0.1,
        uU1=0, sigmaU1=1, uU2=0, sigmaU2=1,
        uV1=0, sigmaV1=1, uV2=0, sigmaV2=1,
        seed=0):
    """ Simulates sigmoid function for microbe-metabolite interations.

    Parameters
    ----------
    num_microbes : int
       Number of microbial species to simulate
    num_metabolites : int
       Number of molecules to simulate
    num_samples : int
       Number of samples to generate
    num_latent_microbes :
       Number of latent microbial dimensions
    num_latent_metabolites
       Number of latent metabolite dimensions
    num_latent_shared
       Number of dimensions in shared representation
    low : float
       Lower bound of gradient
    high : float
       Upper bound of gradient
    microbe_total : int
       Total number of microbial species
    metabolite_total : int
       Total number of metabolite species
    uB : float
       Mean of regression coefficient distribution
    sigmaB : float
       Standard deviation of regression coefficient distribution
    sigmaQ : float
       Standard deviation of error distribution
    uU1 : float
       Mean of microbial input projection coefficient distribution
    sigmaU1 : float
       Standard deviation of microbial input projection
       coefficient distribution
    uU2 : float
       Mean of microbe output projection coefficient distribution
    sigmaU2 : float
       Standard deviation of microbe output projection
       coefficient distribution
    uV1 : float
       Mean of metabolite input projection coefficient distribution
    sigmaU1 : float
       Standard deviation of metabolite input projection
       coefficient distribution
    uV2 : float
       Mean of metabolite output projection coefficient distribution
    sigmaU2 : float
       Standard deviation of metabolite output projection
       coefficient distribution
    seed : float
       Random seed
    Returns
    -------
    microbe_counts : pd.DataFrame
       Count table of microbial counts
    metabolite_counts : pd.DataFrame
       Count table of metabolite counts
    """
    k = num_latent_shared
    state = check_random_state(seed)
    # only have two coefficients
    beta = state.normal(uB, sigmaB, size=(2, k))

    X = np.vstack((np.ones(num_samples),
                   np.linspace(low, high, num_samples))).T

    Q = np.tanh(state.normal(X @ beta, sigmaQ))

    U1 = state.normal(
        uU1, sigmaU1, size=(num_latent_microbes, num_microbes))
    U2 = state.normal(
        uU2, sigmaU2, size=(k, num_latent_microbes))
    V1 = state.normal(
        uV1, sigmaV1, size=(num_latent_metabolites, num_metabolites))
    V2 = state.normal(
        uV2, sigmaV2, size=(k, num_latent_metabolites))

    def multinomial(n, p):
        return np.vstack([np.random.multinomial(n, p[i, :])
                          for i in range(p.shape[0])]).T

    microbe_counts = multinomial(microbe_total, softmax((Q @ U2 @ U1).T))
    metabolite_counts = multinomial(metabolite_total, softmax((Q @ V2 @ V1).T))
    otu_ids = ['OTU_%d' % d for d in range(microbe_counts.shape[1])]
    ms_ids = ['metabolite_%d' % d for d in range(metabolite_counts.shape[1])]
    sample_ids = ['sample_%d' % d for d in range(metabolite_counts.shape[0])]

    microbe_counts = pd.DataFrame(
        microbe_counts, index=sample_ids, columns=otu_ids)
    metabolite_counts = pd.DataFrame(
        metabolite_counts, index=sample_ids, columns=ms_ids)

    return microbe_counts, metabolite_counts, X, Q, U1, U2, V1, V2


def random_multimodal(num_microbes=20, num_metabolites=100, num_samples=100,
                      latent_dim=3, low=-1, high=1,
                      microbe_total=10, metabolite_total=100,
                      uB=0, sigmaB=2, sigmaQ=0.1,
                      uU=0, sigmaU=1, uV=0, sigmaV=1,
                      seed=0):
    """
    Parameters
    ----------
    num_microbes : int
       Number of microbial species to simulate
    num_metabolites : int
       Number of molecules to simulate
    num_samples : int
       Number of samples to generate
    latent_dim :
       Number of latent dimensions
    low : float
       Lower bound of gradient
    high : float
       Upper bound of gradient
    microbe_total : int
       Total number of microbial species
    metabolite_total : int
       Total number of metabolite species
    uB : float
       Mean of regression coefficient distribution
    sigmaB : float
       Standard deviation of regression coefficient distribution
    sigmaQ : float
       Standard deviation of error distribution
    uU : float
       Mean of microbial input projection coefficient distribution
    sigmaU : float
       Standard deviation of microbial input projection
       coefficient distribution
    uV : float
       Mean of metabolite output projection coefficient distribution
    sigmaV : float
       Standard deviation of metabolite output projection
       coefficient distribution
    seed : float
       Random seed

    Returns
    -------
    microbe_counts : pd.DataFrame
       Count table of microbial counts
    metabolite_counts : pd.DataFrame
       Count table of metabolite counts
    """
    state = check_random_state(seed)
    # only have two coefficients
    beta = state.normal(uB, sigmaB, size=(2, num_microbes))
    X = np.vstack((np.ones(num_samples),
                   np.linspace(low, high, num_samples))).T

    microbes = softmax(state.normal(X @ beta, sigmaQ))

    #microbes = softmax(
    #    state.normal(loc=0, scale=sigmaQ,
    #                 size=(num_samples, num_microbes)
    #                )
    #)

    microbes = ilr_inv(state.multivariate_normal(
        mean=np.zeros(num_microbes-1), cov=np.diag([sigmaQ]*(num_microbes-1)),
        size=num_samples)
    )
    Umain = state.normal(
        uU, sigmaU, size=(num_microbes, latent_dim))
    Vmain = state.normal(
        uV, sigmaV, size=(latent_dim, num_metabolites-1))

    Ubias = state.normal(
        uU, sigmaU, size=(num_microbes, 1))
    Vbias = state.normal(
        uV, sigmaV, size=(1, num_metabolites-1))

    U_ = np.hstack(
        (np.ones((num_microbes, 1)), Ubias, Umain))
    V_ = np.vstack(
        (Vbias, np.ones((1, num_metabolites-1)), Vmain))

    phi = np.hstack((np.zeros((num_microbes, 1)), U_ @ V_))
    probs = softmax(phi)
    microbe_counts = np.zeros((num_samples, num_microbes))
    metabolite_counts = np.zeros((num_samples, num_metabolites))
    n1 = microbe_total
    n2 = metabolite_total // microbe_total
    for n in range(num_samples):
        otu = np.random.multinomial(n1, microbes[n, :])
        for i in range(num_microbes):
            ms = np.random.multinomial(otu[i] * n2, probs[i, :])
            metabolite_counts[n, :] += ms
        microbe_counts[n, :] += otu

    otu_ids = ['OTU_%d' % d for d in range(microbe_counts.shape[1])]
    ms_ids = ['metabolite_%d' % d for d in range(metabolite_counts.shape[1])]
    sample_ids = ['sample_%d' % d for d in range(metabolite_counts.shape[0])]

    microbe_counts = pd.DataFrame(
        microbe_counts, index=sample_ids, columns=otu_ids)
    metabolite_counts = pd.DataFrame(
        metabolite_counts, index=sample_ids, columns=ms_ids)

    return microbe_counts, metabolite_counts, X, beta, U_, V_


def ground_truth_edges(microbe_df, metabolite_df):
    """ Hard coded example of edges. """
    interactions = {('theta_f', 'SG'): 1,
                    ('theta_f', 'F'): 1,
                    ('theta_f', 'I'): -1,
                    ('theta_p', 'SA'): 1,
                    ('theta_p', 'P'): 1}
    strains = list(map(lambda x: '_'.join(x.split('_')[:-1]), microbe_df.columns))
    chemicals = list(map(lambda x: '_'.join(x.split('_')[:-1]), metabolite_df.columns))
    edges = []
    for i, otu in enumerate(strains):
        for j, ms in enumerate(chemicals):
            if (otu, ms) not in interactions.keys():
                edges.append((microbe_df.columns[i], metabolite_df.columns[j], 0))
            else:
                direction = interactions[(otu, ms)]
                edges.append((microbe_df.columns[i], metabolite_df.columns[j], direction))
    edges = pd.DataFrame(edges, columns=['microbe', 'metabolite', 'direction'])
    return edges


def random_biofilm(table, uU, sigmaU, uV, sigmaV, sigmaQ,
                   num_microbes, num_metabolites, latent_dim,
                   microbe_total, microbe_kappa,
                   microbe_tau,
                   metabolite_total, metabolite_kappa,
                   metabolite_tau,
                   seed=0):
    """ Generate random biofilm simulation.

    Parameters
    ----------
    table : pd.DataFrame
       Simulated absolute abundances for microbes and metabolites
       along with time and spatial coordinates.
    uU : float
       Mean of microbial input projection coefficient distribution
    sigmaU : float
       Standard deviation of microbial input projection
       coefficient distribution
    uV : float
       Mean of metabolite output projection coefficient distribution
    sigmaV : float
       Standard deviation of metabolite output projection
       coefficient distribution
    sigmaQ : float
       Standard deviation of error distribution
    num_microbes : int
        Number of strains to be represented per microbe
    num_metabolites : int
        Number of chemicals to be represented per metabolite
    microbe_total : int
        Mean total abundance per microbe sample
    microbe_kappa : float
        Dispersion factor for microbes
    microbe_tau : float
        Dispersion factor for total microbial load
    metabolite_total : float
        Mean total intensity per metabolite sample
    metabolite_kappa : float
        Dispersion factor for individual metabolites
    metabolite_tau : float
        Dispersion factor for total metabolite load
    seed : float
       Random seed

    Returns
    -------
    edges : list of tuple
        Edge list of microbes and metabolites specifying the direction of
        interaction.  For instance (x, y, -1) would specify a negative
        correlation between microbe x, metabolite y.  (x, y, 1) would
        specify a positive correlation between microbe x and metabolite y.
    abs_microbes_df : pd.DataFrame
        Relative microbial counts.
    abs_metabolites_df : pd.DataFrame
        Relative metabolite intensities.
    rel_microbes_df : pd.DataFrame
        Relative microbial counts.
    rel_metabolites_df : pd.DataFrame
        Relative metabolite intensities.
    """
    state = check_random_state(seed)
    # Note : this is hard coded
    pairs = [('theta_p', ['P', 'SA']),
             ('theta_f', ['F', 'SG', 'I'])]

    odfs = []
    mdfs = []
    edges = []
    for otu, spectra in pairs:

        # partition the microbes
        microbes_out = partition_microbes(
            num_microbes, sigmaQ, table[otu].values, state
        )

        # partition the metabolites
        for ms in spectra:
            U, V, metabolites_out = partition_metabolites(
                uU, sigmaU, uV, sigmaV,
                num_metabolites,
                latent_dim,
                microbes_out,
                table[ms].values,
                state
            )
            metabolites_df = pd.DataFrame(
                metabolites_out,
                columns=['%s_%d' % (ms, i)
                         for i in range(metabolites_out.shape[1])]
            )
            mdfs.append(metabolites_df)

        microbes_df = pd.DataFrame(
            microbes_out,
            columns=['%s_%d' % (otu, i)
                     for i in range(microbes_out.shape[1])]
        )
        odfs.append(microbes_df)

    ## Helper functions
    # Convert microbial abundances to counts
    def to_counts_f(x):
        n = state.lognormal(np.log(microbe_total), microbe_tau)
        p = x / x.sum()
        return state.poisson(state.lognormal(np.log(n*p), microbe_kappa))

    # Convert metabolite abundances to intensities
    def to_intensities_f(x):
        n = state.lognormal(np.log(metabolite_total), metabolite_tau)
        p = x / x.sum()
        y = state.lognormal(np.log(n*p), metabolite_kappa)
        return y

    abs_microbes_df = pd.concat(odfs, axis=1)
    abs_metabolites_df = pd.concat(mdfs, axis=1)
    abs_microbes_df.index = table.index
    abs_metabolites_df.index = table.index

    rel_microbes_df = abs_microbes_df.apply(to_counts_f, axis=1)
    rel_metabolites_df = abs_metabolites_df.apply(to_intensities_f, axis=1)

    rel_microbes_df = pd.DataFrame(np.vstack(rel_microbes_df.values),
                                   index=abs_microbes_df.index,
                                   columns=abs_microbes_df.columns)

    rel_metabolites_df = pd.DataFrame(np.vstack(rel_metabolites_df.values),
                                      index=abs_metabolites_df.index,
                                      columns=abs_metabolites_df.columns)

    # ground truth edges
    edges = ground_truth_edges(abs_microbes_df, abs_metabolites_df)

    return (edges, abs_microbes_df, abs_metabolites_df,
            rel_microbes_df, rel_metabolites_df)
