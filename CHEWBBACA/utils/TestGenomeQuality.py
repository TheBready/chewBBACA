#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHORS

    Mickael Silva
    email: @mickaelsilva

    Pedro Cerqueira
    github: @pedrorvc

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

"""


import os
import time
from copy import deepcopy
from collections import Counter

import plotly
import numpy as np
import plotly.graph_objs as go


def presence_absence_matrix(matrix):
    """ Converts a matrix from the AlleleCall process
        into a Presence(1)/Absence(0) matrix.

        Args:
            matrix (numpy.ndarray): array with the profiles
            determined by the AlleleCall process of chewBBACCA.
        Returns:
            binary_matrix (numpy.ndarray): array with valid allele
            identifiers (any allele integer identifier, including
            inferred) converted to 1 and missing values (LNF, ASM,
            ALM, PLOT3, PLOT5...) converted to 0.
    """

    binary_matrix = deepcopy(matrix)

    row = 1
    total_rows = binary_matrix.shape[0]
    total_cols = binary_matrix.shape[1]
    while row < total_rows:
        column = 1
        while column < total_cols:
            current = binary_matrix[row, column]
            # convert to integer and change to
            # presence(1)/absence(0)
            try:
                binary_matrix[row, column] = \
                    1 if int(current) > 0 else 0
            # if it cannot be converted to an integer
            # it is either an inferred allele or
            # missing data
            except ValueError:
                binary_matrix[row, column] = \
                    1 if 'INF-' in current else 0

            column += 1
        row += 1

    return binary_matrix


def presence_absence_iterator(matrix, threshold, vector):
    """ Determines the number of genes present in 0%, 95%, 99%,
        99.5% and 100% of genomes, excluding the genomes that
        have a number of missing genes greater than a defined
        threshold.

        Args:
            matrix (numpy.ndarray): array with the profiles
            determined by the AlleleCall process of chewBBACCA
            in a binary format (Presence/Absence matrix).
            threshold (int): maximum acceptable number of
            missing genes. Genomes with a number of missing
            genes above this threshold will be excluded from
            the calculations to determine the number of present
            genes.
            vector (list): list with 7 sublists. The sublists
            will store (in order): number of genomes included
            in the calculation of present genes (respect the
            missing genes threshold), number of genes present
            in 95% of the genomes, number of genes present in
            99% of the genomes, number of genes present in 
            99.5% of the genomes, number of genes present in 
            0% of the genomes, number of genomes excluded from
            the analysis and number of genes present in 100% 
            of the genomes.

        Returns:
            exclude_genomes (list): genomes that had a number
            of missing loci above the threshold.
            vector (list): updated vector with the results for
            the current threshold.
            present_genes (list): identifiers of genes that are
            present in at least 95% of the genomes.
    """

    iteration_matrix = deepcopy(matrix)

    genomes = iteration_matrix[1:, :1]

    # this list includes 'FILE' header from matrix
    genes = iteration_matrix[0]

    column = 1
    totals = []
    # stores list of genomes that had a number of missing loci
    # greater than the threshold
    exclude_genomes = []
    # stores genes that were not in any of the genomes
    # stores genes that were in all genomes
    genes_0 = []
    genes_95 = 0
    genes_99 = 0
    genes_995 = 0
    genes_100 = []
    present_genes = []

    # iterate over all loci columns
    total_rows = iteration_matrix.shape[0]
    total_cols = iteration_matrix.shape[1]
    while column < total_cols:
        # start at row of first genome
        # row 0 is 'FILE'
        row = 1
        missing_genes = 0
        genome_miss = []
        # iterate over genomes rows for each locus column
        while row < total_rows:
            current_cell = int(iteration_matrix[row, column])
            # if cell value is '0', the gene was not found in
            # that genome
            if current_cell == 0:
                # increment number of times gene is not found
                missing_genes += 1
                # store genome index as bad genome (does not
                # have the gene)
                genome_miss.append(genomes[int(row - 1)][0])

            # increment to check next genome
            row += 1

        # after iterating over all rows/genomes for a locus column
        # calculate percentage of genomes that have the locus
        gene_presence = float((total_rows-1) - missing_genes) / float(total_rows-1)

        # stores values to show in plotly plot
        if gene_presence >= 0.95:
            genes_95 += 1
        if gene_presence >= 0.99:
            genes_99 += 1
        if gene_presence >= 0.995:
            genes_995 += 1

        if len(genomes) > 500:
            xthreshold = 0.95
        elif len(genomes) < 20:
            xthreshold = 0.90
        # most of the times, it will be 0.95
        else:
            xthreshold = 0.95

        # selects genes present in at least 95% of genomes
        if gene_presence >= xthreshold:
            # list that stores values to show in 'GENES_95%.txt'
            present_genes.append(genes[column])
            if int(gene_presence) == 1:
                genes_100.append(genes[column])
        elif gene_presence == 0:
            genes_0.append(genes[column])

        if gene_presence > xthreshold and len(genome_miss) > 0:
            for genome in genome_miss:
                # get identifiers of genomes that did not have the
                # gene/locus
                exclude_genomes.append(genome)

        # we need to subtract 1 to total_rows shape because 'FILE'
        # header is being included in genome list
        totals.append(float(float((total_rows-1) - missing_genes) / float(total_rows-1)))
        column += 1

    counter = Counter(exclude_genomes)
    most_common = counter.most_common()

    exclude_genomes = [genome[0] for genome in most_common
                       if int(genome[1]) > threshold]

    # number of used genomes
    vector[0].append(len(genomes))

    # number of loci at 95%
    vector[1].append(genes_95)

    # number of loci at 99%
    vector[2].append(genes_99)

    # number of loci at 99.5%
    vector[3].append(genes_995)

    # number of loci at 0%
    vector[4].append(len(genes_0))

    # number of genomes to be removed
    # if the sublist already has a value, sum the new value
    try:
        vector[5].append(len(exclude_genomes) + ((vector[5])[-1]))
    # if it is empty just append
    except Exception:
        vector[5].append(len(exclude_genomes))

    # number of loci at 100%
    vector[6].append(len(genes_100))

    return [exclude_genomes, vector, present_genes]


def remove_genomes(matrix, rm_genomes):
    """ Removes matrix rows corresponding to
        genomes listed in the 'rm_genomes' argument.

        Args:
            matrix (numpy.ndarray): matrix with the
            profiles determined by the AlleleCall
            process (accepts original or presence/
            absence matrix).
            rm_genomes (list): identifiers of genomes
            whose rows will be removed.

        Returns:
            prunned_matrix (numpy.ndarray): input matrix
            without the rows corresponding to genomes
            listed in the 'rm_genomes' argument.
    """

    prunned_matrix = deepcopy(matrix)

    genomes = (prunned_matrix[1:, :1]).tolist()

    i = 0
    for genome in genomes:
        i += 1
        if genome[0] in rm_genomes:
            prunned_matrix = np.delete(prunned_matrix, i, axis=0)
            i -= 1

    return prunned_matrix


def threshold_info(matrix, iterations, threshold):
    """ Applies an iterative approach to determine stable
        values for the number of genes present in 0%, 95%, 99%,
        99.5% and 100% of genomes, excluding the genomes that
        have a number of missing genes greater than a defined
        threshold.

        Args:
            matrix (numpy.ndarray): matrix with the profiles
            determined by the AlleleCall process in the format
            of a presence/absence matrix.
            iterations (int): maximum number of iterations that
            will be used to reach stable values of genomes used
            and present genes.
            threshold (int): maximum acceptable number of missing
            genes. Genomes with a number of missing genes greater
            that this threshold will be excluded at each iteration.

        Returns:
            stats_vector (list): list with 7 sublists. The sublists
            will store (in order): number of genomes included
            in the calculation of present genes (respect the
            missing genes threshold), number of genes present
            in 95% of the genomes, number of genes present in
            99% of the genomes, number of genes present in
            99.5% of the genomes, number of genes present in
            0% of the genomes, number of genomes excluded from
            the analysis and number of genes present in 100%
            of the genomes.
            stable_iteration (int): the iteration at which the
            number of used genomes and present genes stabilized.
            removed_genomes (list): genomes excluded from the
            analysis at the current threshold.
            iterations_genes (list): identifiers of genes that
            are present in 95% of the genomes.
    """

    i = 0
    removed_genomes = []
    genomes_to_remove = []

    stats_vector = [[] for x in range(7)]

    stable = False
    present_genes = []
    iterations_genes = []
    stable_iteration = None
    removed_genomes_count = 0
    while i <= iterations or not stable:

        # genomes to remove because they had
        # a number of missing loci greater than threshold
        for g in genomes_to_remove:
            removed_genomes.append(g)

        # keep iterating if the number of genomes to remove keeps
        # increasing
        if len(removed_genomes) > removed_genomes_count and i > 0:
            removed_genomes_count = len(removed_genomes)
        elif stable_iteration is None and i > 0:
            stable = True
            stable_iteration = i
            iterations_genes.append(present_genes)
        if not stable:
            # remove genomes rows from matrix and re-calculate
            # the gene presence percentage
            if len(genomes_to_remove) > 0:
                matrix = remove_genomes(matrix, genomes_to_remove)

            genomes_to_remove, stats_vector, present_genes = \
                presence_absence_iterator(matrix,
                                          threshold,
                                          stats_vector)

        # increment iterations
        i += 1

    return [stats_vector, stable_iteration,
            removed_genomes, iterations_genes]


def scatter_tracer(x_data, y_data, tracer_name, tracer_mode,
                   ref_yaxis, marker_symbol, marker_size,
                   marker_color, line_dash):
    """ Creates a tracer object for a scatter plot.

        Args:
            x_data (list): xaxis values corresponding
            to the maximum number of missing genes
            threshold values.
            y_data (list): yaxis values corresponding
            to either the number of genomes included
            in the analysis at each threshold or to
            the number of genes present in 95%, 99%,
            99.5%, or 100% of the genomes at each
            threshold.
            tracer_name (str): name to show in the
            plot legend.
            tracer_mode (str): type of symbols used
            used to represent data (lines, markers,
            lines+markers...).
            ref_yaxis (str): the yaxis that will be
            used as reference.
            marker_symbol (str): the type of symbol
            used for the markers.
            marker_size (int): the size of the marker
            symbol.
            marker_color (str): color of the markers
            and of the line (if any).
            line_dash (str): type of line (solid,
            dash...).

        Returns:
            tracer (plotly.graph_objs.Scatter): a
            Plotly tracer with the data for a scatter
            plot (a group of data points).
    """

    tracer = go.Scatter(x=x_data,
                        y=y_data,
                        name=tracer_name,
                        mode=tracer_mode,
                        yaxis=ref_yaxis,
                        marker=dict(symbol=marker_symbol,
                                    size=marker_size,
                                    color=marker_color),
                        line=dict(dash=line_dash)
                        )

    return tracer


def main(input_path, num_iterations,
         missing_loci_threshold, step, out_folder):

    # create output directory, if it does not exist
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    start_time = '\nStarting Script at: {0}'.format(time.strftime('%H:%M:%S-%d/%m/%Y'))

    # define starting threshold
    # starting with maximum number of missing loci of 0
    threshold = 0
    results = []
    stable_iterations = []
    threshold_list = list(range(0, missing_loci_threshold+step, step))

    # store the list of genomes removed per threshold
    removed_genomes_file = os.path.join(out_folder, 'removedGenomes.txt')
    with open(removed_genomes_file, 'w') as rgf:
        rgf.write('Threshold\tRemoved_genomes\n')

    # store the list of schema genes/loci per threshold
    present_genes_file = os.path.join(out_folder, 'Genes_95%.txt')
    with open(present_genes_file, 'w') as pgf:
        pgf.write('Threshold\tPresent_genes\n')

    print('\nImporting matrix and converting to Presence/Absence matrix...')
    # read AlleleCall matrix
    original_matrix = np.genfromtxt(input_path, delimiter='\t',
                                    dtype=None, encoding=None)

    # copy and work with copied version
    copy_matrix = np.copy(original_matrix)
    # create presence/absence matrix for easier processing
    copy_matrix = presence_absence_matrix(copy_matrix)
    print('Done.')

    # determine genomes and genes for every threshold
    print('\nMissing genes thresholds:'
          '\nMinimun: {0}'
          '\nMaximum: {1}'
          '\nStep: {2}'.format(threshold, missing_loci_threshold, step))

    print('\nDetermining number of genes present in '
          '95%, 99%, 99.5% and 100% of genomes for all threshold values...\n')

    print('-'*67)
    print('{0:^9}  {1:^9}  {2:^7}  '
          '{3:^7}  {4:^7}  {5:^8}  {6:^8}'.format('THRESHOLD',
                                                  'ITERATION',
                                                  'GENOMES',
                                                  'GENES95',
                                                  'GENES99',
                                                  'GENES995',
                                                  'GENES100'))
    print('-'*67)

    total_genomes = original_matrix.shape[0] - 1
    while threshold <= missing_loci_threshold:

        result, stable_iteration, removed_genomes, iterations_genes = threshold_info(copy_matrix,
                                                                                     num_iterations,
                                                                                     threshold)

        stable_iterations.append(stable_iteration)

        # print threshold line
        print('{0:^9}  {1:^9}  {2:^7}  '
              '{3:^7}  {4:^7}  {5:^8}  {6:^8}'.format(threshold,
                                                      stable_iteration,
                                                      total_genomes - len(removed_genomes),
                                                      result[1][-1],
                                                      result[2][-1],
                                                      result[3][-1],
                                                      result[6][-1]))

        results.append(result)

        # save results to files
        with open(removed_genomes_file, 'a') as rgf:
            rm_genomes_text = '{0}\t{1}\n'.format(threshold,
                                                  ' '.join(map(str, removed_genomes)))
            rgf.write(rm_genomes_text)

        with open(present_genes_file, 'a') as pgf:
            pgf.write('{0}\t'.format(threshold))

            for x in iterations_genes:
                pgf.write((' '.join(map(str, x))) + '\n')

        # increment maximum number of missing loci
        threshold += step

    print('-'*67)

    # plot legend labels
    labels = ['Genomes used',
              'Loci present in 95% genomes',
              'Loci present in 99% genomes',
              'Loci present in 99.5% genomes',
              'Loci present in 100% genomes']

    colors = ['#3690c0', '#ec7014', '#66bd63', '#807dba', '#fdbb84']
    # xaxis data is the list of thresholds
    x_data = threshold_list

    # number of genomes used per threshold
    y_genomes = [res[0][-1] for res in results]
    # number of genes per threshold and per
    # genome presence percentage
    y_95 = [res[1][-1] for res in results]
    y_99 = [res[2][-1] for res in results]
    y_995 = [res[3][-1] for res in results]
    y_100 = [res[6][-1] for res in results]

    # group all yaxis datasets into list
    y_datasets = [y_genomes, y_95, y_99, y_995, y_100]

    # create all tracers
    tracers = []
    for d in range(len(y_datasets)):
        # tracer with used genomes data
        if d == 0:
            tracer = scatter_tracer(x_data,
                                    y_datasets[d],
                                    labels[d],
                                    'lines+markers',
                                    'y2',
                                    'diamond-dot',
                                    10,
                                    colors[d],
                                    'solid')
        # tracers for number of genes per threshold
        else:
            tracer = scatter_tracer(x_data,
                                    y_datasets[d],
                                    labels[d],
                                    'lines+markers',
                                    'y1',
                                    'circle',
                                    10,
                                    colors[d],
                                    'dash')

        tracers.append(tracer)

    # define layout attributes
    fig_layout = go.Layout(title='Test genomes quality',
                           xaxis=dict(title='Threshold',
                                      showline=True,
                                      mirror=True,
                                      linecolor='#EBEBEB',
                                      ticks='outside',
                                      tickcolor='#EBEBEB',
                                      showgrid=True,
                                      gridcolor='rgb(255,255,255)',
                                      range=[-5, missing_loci_threshold+5]),
                           yaxis=dict(title='Number of loci',
                                      showline=True,
                                      linecolor='#EBEBEB',
                                      ticks='outside',
                                      tickcolor='#EBEBEB',
                                      showgrid=True,
                                      gridcolor='rgb(255,255,255)'),
                           yaxis2=dict(title='Number of genomes',
                                       showline=True,
                                       linecolor='#EBEBEB',
                                       ticks='outside',
                                       tickcolor='#EBEBEB',
                                       showgrid=True,
                                       gridcolor='rgb(255,255,255)',
                                       overlaying='y',
                                       side='right'),
                           plot_bgcolor='#EBEBEB'
                           )

    fig = go.Figure(data=tracers, layout=fig_layout)
    plot_file = os.path.join(out_folder, 'GenomeQualityPlot.html')
    plotly.offline.plot(fig, filename=plot_file)

    print(start_time)
    print('Finished Script at: {0}'.format(time.strftime('%H:%M:%S-%d/%m/%Y')))


if __name__ == "__main__":

    main()
