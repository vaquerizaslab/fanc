import argparse
import textwrap
import os


def klot_parser():
    parser = argparse.ArgumentParser(
        description="klot plotting tool for kaic",
        usage='''
            klot [<klot global parameters>] <region> [<region> ...]
                  --plot <plot data file(s)> [<plot parameters>] [...]

            Run klot --plot -t <plot type> -h for help on a specific subplot.
            '''
    )

    parser.add_argument(
        'regions',
        nargs='*',
        default=[],
        help='List of region selectors (<chr>:<start>-<end>) or files with region information (BED, GTF, ...).'
    )

    parser.add_argument(
        '-o', '--output', dest='output',
        help='''Suppresses interactive plotting window and redirects plot to file.
                Specify path to file when plotting a single region,
                and path to a folder for plotting multiple regions.'''
    )

    parser.add_argument(
        '-s', '--script', dest='script',
        help='''Use a script file to define plot.'''
    )

    parser.add_argument(
        '-p', '--plot', dest='plot',
        action='append',
        help='''New plot, type will be chosen automatically by file type, unless '-t' is provided.'''
    )

    parser.add_argument(
        '-n', '--name', dest='name',
        default='',
        help='''Plot name to be used as prefix when plotting multiple regions. Is ignored for single region and
                interactive plot.'''
    )

    parser.add_argument(
        '--height', dest='height',
        type=int,
        help='''Height of the figure in inches. Default is proportional to figure width, dependent on number
                of subplots.'''
    )

    parser.add_argument(
        '--width', dest='width',
        type=int,
        default=6,
        help='''Width of the figure in inches. Default: 6'''
    )

    parser.add_argument(
        '-w', '--window-size', dest='window_size',
        type=int,
        help='''Plotting region size in base pairs. If provided, the actual size of the given region is
                ignored and instead a region
                <chromosome>:<region center - window size/2>-<region center + window size/2> will be plotted.'''
    )

    parser.add_argument(
        '-vs', '--vertical-space', dest='hspace',
        type=float,
        default=.5,
        help='''Vertical distance between plots in fraction of figure.'''
    )

    parser.add_argument(
        '--invert-x', dest='invert_x',
        action='store_true',
        help='''Invert x-axis for this plot'''
    )
    parser.set_defaults(invert_x=False)

    parser.add_argument(
        '-v', '--version', dest='print_version',
        action='store_true',
        help='''Print version information'''
    )
    parser.set_defaults(print_version=False)
    return parser


def type_parser():
    parser = argparse.ArgumentParser(
        description="klot subplot identifier",
        add_help=False,
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        'data',
        nargs='*',
        help='Data to be plotted in subplot.'
    )

    parser.add_argument(
        '-t', '--type', dest='type',
        help=textwrap.dedent('''\
            Manually specify subplot type. Options:
            hic          Hi-C plot, cropped triangle style
            hic2d        Hi-C plot, matrix style
            hicsplit     Hi-C vs Hi-C plot, split matrix
            hicvhic      Hi-C vs Hi-C plot, matrices "mirrored"
            fc           Fold-change plot, cropped triangle style
            hicvfc       Hi-C vs fold-change plot, matrices "mirrored"
            array        Array "flame" plot (e.g. insulation index)
            region       Bar plot with region score (e.g. BED)
            line         Line plot with values per region
            bigwig       Plot BigWig files
            gene         Genes plot (exon/intron structure)
        ''')
    )
    return parser


def subplot_parser():
    parser = argparse.ArgumentParser(
        description="klot subplot identifier",
    )

    parser.add_argument(
        '--aspect-ratio', dest='aspect_ratio',
        type=float,
        help='''Aspect ratio of this panel. Default is determined by figure type (usually 1.0).'''
    )

    parser.add_argument(
        '--title', dest='title',
        default='',
        help='''Title of this plot.'''
    )

    parser.add_argument(
        '--fix-chromosome', dest='fix_chromosome',
        action='store_true',
        help='''Fix chromosome identifier for this plot (add or remove 'chr' as required)'''
    )
    parser.set_defaults(fix_chromosome=False)

    parser.add_argument(
        '--hide-x', dest='hide_x',
        action='store_true',
        help='''Hide x-axis for this plot'''
    )
    parser.set_defaults(hide_x=False)

    return parser


def fc(parameters):
    parser = subplot_parser()
    parser.description = '''Fold-change plot.'''

    parser.add_argument(
        'fc_matrix',
        help='''Fold-change matrix.'''
    )

    parser.add_argument(
        '-vmin', '--minimum-value', dest='vmin',
        type=float,
        help='''Minimum value assigned the first color in the colorbar.'''
    )

    parser.add_argument(
        '-vmax', '--maximum-value', dest='vmax',
        type=float,
        help='''Maximum value assigned the last color in the colorbar.'''
    )

    parser.add_argument(
        '-d', '--maximum-distance', dest='max_dist',
        type=int,
        help='''Maximum distance between two points after which triangle will be truncated.'''
    )

    parser.add_argument(
        '-l', '--log', dest='log',
        action='store_true',
        help='''Log-transform heatmap values'''
    )
    parser.set_defaults(log=False)

    parser.add_argument(
        '-S', '--no-symmetry', dest='symmetry',
        action='store_false',
        help='''Do not plot colormap symmetrical around 0.'''
    )
    parser.set_defaults(symmetry=True)

    parser.add_argument(
        '-r', '--range-slider', dest='adjust_range',
        action='store_true',
        help='''Add vmax/vmin slider to plot'''
    )
    parser.set_defaults(adjust_range=False)

    parser.add_argument(
        '-c', '--colormap', dest='colormap',
        default='RdBu_r',
        help='''Matplotlib colormap (default: RdBu_r)'''
    )

    parser.add_argument(
        '-C', '--no-colorbar', dest='show_colorbar',
        action='store_false',
        help='''Do not show colorbar in plot'''
    )
    parser.set_defaults(show_colorbar=True)

    args = parser.parse_args(parameters)

    import kaic
    import kaic.plotting as kplt
    matrix = kaic.FoldChangeMatrix(os.path.expanduser(args.fc_matrix), mode='r')
    norm = "lin" if not args.log else "log"
    colorbar_symmetry = 0 if args.symmetry else None
    return kplt.HicPlot(matrix, colormap=args.colormap, max_dist=args.max_dist, norm=norm, vmin=args.vmin,
                        vmax=args.vmax, show_colorbar=args.show_colorbar, adjust_range=args.adjust_range,
                        colorbar_symmetry=colorbar_symmetry), args


def hic_parser():
    parser = subplot_parser()
    parser.description = '''Hi-C plot.'''

    parser.add_argument(
        'hic',
        help='''Hi-C object.'''
    )

    parser.add_argument(
        '-vmin', '--minimum-value', dest='vmin',
        type=float,
        help='''Minimum value assigned the first color in the colorbar.'''
    )

    parser.add_argument(
        '-vmax', '--maximum-value', dest='vmax',
        type=float,
        help='''Maximum value assigned the last color in the colorbar.'''
    )

    parser.add_argument(
        '-d', '--maximum-distance', dest='max_dist',
        type=int,
        help='''Maximum distance between two points after which triangle will be truncated.'''
    )

    parser.add_argument(
        '-l', '--log', dest='log',
        action='store_true',
        help='''Log-transform heatmap values'''
    )
    parser.set_defaults(log=False)

    parser.add_argument(
        '-r', '--range-slider', dest='adjust_range',
        action='store_true',
        help='''Add vmax/vmin slider to plot'''
    )
    parser.set_defaults(adjust_range=False)

    parser.add_argument(
        '-c', '--colormap', dest='colormap',
        help='''Matplotlib colormap'''
    )

    parser.add_argument(
        '-C', '--no-colorbar', dest='show_colorbar',
        action='store_false',
        help='''Do not show colorbar in plot'''
    )
    parser.set_defaults(show_colorbar=True)
    return parser


def hic(parameters):
    parser = hic_parser()
    args = parser.parse_args(parameters)

    import kaic
    from kaic.config import config
    import kaic.plotting as kplt
    colormap = config.colormap_hic if args.colormap is None else args.colormap

    matrix = kaic.load_hic(os.path.expanduser(args.hic), mode='r')

    norm = "lin" if not args.log else "log"
    return kplt.HicPlot(matrix, colormap=colormap, max_dist=args.max_dist, norm=norm, vmin=args.vmin,
                        vmax=args.vmax, show_colorbar=args.show_colorbar, adjust_range=args.adjust_range), args


def hic2d_parser():
    parser = subplot_parser()
    parser.description = '''Hi-C 2D plot.'''

    parser.add_argument(
        'hic',
        help='''Hi-C object.'''
    )

    parser.add_argument(
        '-vmin', '--minimum-value', dest='vmin',
        type=float,
        help='''Minimum value assigned the first color in the colorbar.'''
    )

    parser.add_argument(
        '-vmax', '--maximum-value', dest='vmax',
        type=float,
        help='''Maximum value assigned the last color in the colorbar.'''
    )

    parser.add_argument(
        '-l', '--log', dest='log',
        action='store_true',
        help='''Log-transform heatmap values'''
    )
    parser.set_defaults(log=False)

    parser.add_argument(
        '-r', '--range-slider', dest='adjust_range',
        action='store_true',
        help='''Add vmax/vmin slider to plot'''
    )
    parser.set_defaults(log=False)

    parser.add_argument(
        '-c', '--colormap', dest='colormap',
        help='''Matplotlib colormap'''
    )

    parser.add_argument(
        '-C', '--no-colorbar', dest='show_colorbar',
        action='store_false',
        help='''Do not show colorbar in plot'''
    )
    parser.set_defaults(show_colorbar=True)
    return parser


def hic2d(parameters):
    parser = hic2d_parser()

    args = parser.parse_args(parameters)
    norm = "lin" if not args.log else "log"

    import kaic
    from kaic.config import config
    import kaic.plotting as kplt
    colormap = config.colormap_hic if args.colormap is None else args.colormap

    matrix = kaic.load_hic(os.path.expanduser(args.hic), mode='r')
    return kplt.HicPlot2D(matrix, colormap=colormap, norm=norm, vmin=args.vmin, vmax=args.vmax,
                          show_colorbar=args.show_colorbar, adjust_range=args.adjust_range), args


def hicsplit_parser():
    parser = subplot_parser()
    parser.description = '''Hi-C vs Hi-C split plot.'''

    parser.add_argument(
        'hic_top',
        help='''Top Hi-C object.'''
    )

    parser.add_argument(
        'hic_bottom',
        help='''Bottom Hi-C object.'''
    )

    parser.add_argument(
        '-vmin', '--minimum-value', dest='vmin',
        type=float,
        help='''Minimum value assigned the first color in the colorbar.'''
    )

    parser.add_argument(
        '-vmax', '--maximum-value', dest='vmax',
        type=float,
        help='''Maximum value assigned the last color in the colorbar.'''
    )

    parser.add_argument(
        '-l', '--log', dest='log',
        action='store_true',
        help='''Log-transform heatmap values'''
    )
    parser.set_defaults(log=False)

    parser.add_argument(
        '-c', '--colormap', dest='colormap',
        help='''Matplotlib colormap'''
    )

    parser.add_argument(
        '-C', '--no-colorbar', dest='show_colorbar',
        action='store_false',
        help='''Do not show colorbar in plot'''
    )
    parser.set_defaults(show_colorbar=True)

    parser.add_argument(
        '-S', '--no-scaling', dest='scaling',
        action='store_false',
        help='''Do not scale matrix values. By default, the matrices are scaled so that they sum up to the same
                    number of contacts.'''
    )
    parser.set_defaults(scaling=True)
    return parser


def hicsplit(parameters):
    parser = hicsplit_parser()

    args = parser.parse_args(parameters)

    import kaic
    from kaic.config import config
    import kaic.plotting as kplt
    colormap = config.colormap_hic if args.colormap is None else args.colormap

    matrix_bottom = kaic.load_hic(os.path.expanduser(args.hic_bottom), mode='r')
    matrix_top = kaic.load_hic(os.path.expanduser(args.hic_top), mode='r')

    norm = "lin" if not args.log else "log"
    sp = kplt.HicComparisonPlot2D(matrix_top, matrix_bottom, colormap=colormap, norm=norm, vmin=args.vmin,
                                  vmax=args.vmax, scale_matrices=args.scaling, show_colorbar=args.show_colorbar)
    return sp, args


def hicvhic_parser():
    parser = subplot_parser()
    parser.description = '''Hi-C vs Hi-C plot.'''

    parser.add_argument(
        'hic_upper',
        help='''Upper Hi-C object.'''
    )

    parser.add_argument(
        'hic_lower',
        help='''Lower Hi-C object.'''
    )

    parser.add_argument(
        '-uvmin', '--minimum-value-upper', dest='vmin_upper',
        type=float,
        help='''Minimum value assigned the first color in the colorbar of upper plot.'''
    )

    parser.add_argument(
        '-uvmax', '--maximum-value-upper', dest='vmax_upper',
        type=float,
        help='''Maximum value assigned the last color in the colorbar of upper plot.'''
    )

    parser.add_argument(
        '-lvmin', '--minimum-value-lower', dest='vmin_lower',
        type=float,
        help='''Minimum value assigned the first color in the colorbar of lower plot.'''
    )

    parser.add_argument(
        '-lvmax', '--maximum-value-lower', dest='vmax_lower',
        type=float,
        help='''Maximum value assigned the last color in the colorbar of lower plot.'''
    )

    parser.add_argument(
        '-d', '--maximum-distance', dest='max_dist',
        type=int,
        help='''Maximum distance between two points after which triangle will be truncated.'''
    )

    parser.add_argument(
        '-ul', '--log-upper', dest='log_upper',
        action='store_true',
        help='''Log-transform heatmap values of upper plot'''
    )
    parser.set_defaults(log_upper=False)

    parser.add_argument(
        '-ll', '--log-lower', dest='log_lower',
        action='store_true',
        help='''Log-transform heatmap values of lower plot'''
    )
    parser.set_defaults(log_lower=False)

    parser.add_argument(
        '-uc', '--colormap-upper', dest='colormap_upper',
        help='''Matplotlib colormap for upper plot'''
    )

    parser.add_argument(
        '-lc', '--colormap-lower', dest='colormap_lower',
        help='''Matplotlib colormap for lower plot'''
    )

    parser.add_argument(
        '-C', '--no-colorbars', dest='show_colorbar',
        action='store_false',
        help='''Do not show colorbars in plot'''
    )
    parser.set_defaults(show_colorbar=True)
    return parser


def hicvhic(parameters):
    parser = hicvhic_parser()

    args = parser.parse_args(parameters)

    import kaic
    from kaic.config import config
    import kaic.plotting as kplt
    colormap_lower = config.colormap_hic if args.colormap_lower is None else args.colormap_lower
    colormap_upper = config.colormap_hic if args.colormap_upper is None else args.colormap_upper

    matrix_upper = kaic.load_hic(os.path.expanduser(args.hic_upper), mode='r')
    matrix_lower = kaic.load_hic(os.path.expanduser(args.hic_lower), mode='r')

    norm_upper = "lin" if not args.log_upper else "log"
    upper_plot = kplt.HicPlot(matrix_upper, colormap=colormap_upper, max_dist=args.max_dist, norm=norm_upper,
                              vmin=args.vmin_upper, vmax=args.vmax_upper, show_colorbar=args.show_colorbar,
                              adjust_range=False)
    norm_lower = "lin" if not args.log_lower else "log"
    lower_plot = kplt.HicPlot(matrix_lower, colormap=colormap_lower, max_dist=args.max_dist,
                              norm=norm_lower, vmin=args.vmin_lower, vmax=args.vmax_lower,
                              show_colorbar=args.show_colorbar, adjust_range=False)
    vsp = kplt.VerticalSplitPlot(upper_plot, lower_plot)
    return vsp, args


def hicvfc_parser():
    parser = subplot_parser()
    parser.description = '''Hi-C vs fold-change plot.'''

    parser.add_argument(
        'hic',
        help='''Hi-C object.'''
    )

    parser.add_argument(
        'fold_change',
        help='''Fold-change object.'''
    )

    parser.add_argument(
        '-vmin', '--minimum-value-hic', dest='vmin_hic',
        type=float,
        help='''Minimum value assigned the first color in the colorbar of Hi-C plot.'''
    )

    parser.add_argument(
        '-vmax', '--maximum-value-hic', dest='vmax_hic',
        type=float,
        help='''Maximum value assigned the last color in the colorbar of Hi-C plot.'''
    )

    parser.add_argument(
        '-fvmin', '--minimum-value-fc', dest='vmin_fc',
        type=float,
        help='''Minimum value assigned the first color in the colorbar of fold-change plot.'''
    )

    parser.add_argument(
        '-fvmax', '--maximum-value-fc', dest='vmax_fc',
        type=float,
        help='''Maximum value assigned the last color in the colorbar of fold-change plot.'''
    )

    parser.add_argument(
        '-d', '--maximum-distance', dest='max_dist',
        type=int,
        help='''Maximum distance between two points after which triangle will be truncated.'''
    )

    parser.add_argument(
        '-hl', '--log-hic', dest='log_hic',
        action='store_true',
        help='''Log-transform heatmap values of Hi-C plot'''
    )
    parser.set_defaults(log_hic=False)

    parser.add_argument(
        '-fl', '--log-fc', dest='log_fc',
        action='store_true',
        help='''Log-transform heatmap values of fold-change plot'''
    )
    parser.set_defaults(log_fc=False)

    parser.add_argument(
        '-hc', '--colormap-hic', dest='colormap_hic',
        help='''Matplotlib colormap for Hi-C plot'''
    )

    parser.add_argument(
        '-fc', '--colormap-fc', dest='colormap_fc',
        default='RdBu_r',
        help='''Matplotlib colormap for fold-change plot (default: RdBu_r)'''
    )

    parser.add_argument(
        '-C', '--no-colorbars', dest='show_colorbar',
        action='store_false',
        help='''Do not show colorbars in plot'''
    )
    parser.set_defaults(show_colorbar=True)

    parser.add_argument(
        '-S', '--no-symmetry', dest='symmetry',
        action='store_false',
        help='''Do not plot colormap symmetrical around 0.'''
    )
    parser.set_defaults(symmetry=True)

    parser.add_argument(
        '-i', '--invert', dest='invert',
        action='store_true',
        help='''Invert plot (fold-change on top)'''
    )
    parser.set_defaults(invert=False)
    return parser


def hicvfc(parameters):
    parser = hicvfc_parser()

    args = parser.parse_args(parameters)

    import kaic
    from kaic.config import config
    import kaic.plotting as kplt

    colormap_hic = config.colormap_hic if args.colormap_hic is None else args.colormap_hic

    hic = kaic.load_hic(os.path.expanduser(args.hic), mode='r')
    fc = kaic.FoldChangeMatrix(os.path.expanduser(args.fold_change), mode='r')

    norm_hic = "lin" if not args.log_hic else "log"
    hic_plot = kplt.HicPlot(hic, colormap=colormap_hic, max_dist=args.max_dist, norm=norm_hic,
                            vmin=args.vmin_hic, vmax=args.vmax_hic, show_colorbar=args.show_colorbar,
                            adjust_range=False)
    norm_fc = "lin" if not args.log_fc else "log"
    colorbar_symmetry = 0 if args.symmetry else None
    fc_plot = kplt.HicPlot(fc, colormap=args.colormap_fc, max_dist=args.max_dist,
                           norm=norm_fc, vmin=args.vmin_fc, vmax=args.vmax_fc,
                           show_colorbar=args.show_colorbar, adjust_range=False,
                           colorbar_symmetry=colorbar_symmetry)
    if args.invert:
        vsp = kplt.VerticalSplitPlot(fc_plot, hic_plot)
    else:
        vsp = kplt.VerticalSplitPlot(hic_plot, fc_plot)

    return vsp, args


def array_parser():
    parser = subplot_parser()
    parser.description = '''Hi-C plot.'''

    parser.add_argument(
        'array',
        help='''Array object, e.g. InsulationIndex, DirectionalityIndex, ... .'''
    )

    parser.add_argument(
        '-vmin', '--minimum-value', dest='vmin',
        type=float,
        help='''Minimum value assigned the first color in the colorbar.'''
    )

    parser.add_argument(
        '-vmax', '--maximum-value', dest='vmax',
        type=float,
        help='''Maximum value assigned the last color in the colorbar.'''
    )

    parser.add_argument(
        '-l', '--log', dest='log',
        action='store_true',
        help='''Log-transform heatmap values'''
    )
    parser.set_defaults(log=False)

    parser.add_argument(
        '-m', '--matrix-height', dest='matrix_height',
        type=int,
        help='''Matrix height in (steps/bins)'''
    )

    parser.add_argument(
        '-r', '--range', dest='range',
        nargs=2,
        type=int,
        help='''Range of y-values to plot (<min> <max> inclusive)'''
    )

    parser.add_argument(
        '-f', '--fields', dest='fields',
        type=str,
        nargs='+',
        help='''List of specific fields to plot.'''
    )

    parser.add_argument(
        '-c', '--colormap', dest='colormap',
        default='RdBu_r',
        help='''Matplotlib colormap (default: RdBu_r)'''
    )

    parser.add_argument(
        '-C', '--no-colorbar', dest='show_colorbar',
        action='store_false',
        help='''Do not show colorbar in plot'''
    )
    parser.set_defaults(show_colorbar=True)

    parser.add_argument(
        '-S', '--no-symmetry', dest='symmetry',
        action='store_false',
        help='''Do not plot colormap symmetrical around 0.'''
    )
    parser.set_defaults(symmetry=True)

    parser.add_argument(
        '-rc', '--replacement-color', dest='replacement_color',
        default='grey',
        help='''Color to replace missing values. Default: grey'''
    )
    return parser


def array(parameters):
    parser = array_parser()

    args = parser.parse_args(parameters)

    import kaic
    import kaic.plotting as kplt

    array = kaic.load(os.path.expanduser(args.array), mode='r')
    norm = "linear" if not args.log else "log"

    data_selection = None
    if args.range is not None:
        data_selection = []
        for i, y in enumerate(array.y_values):
            if args.range[0] <= y <= args.range[1]:
                data_selection.append(i)
    elif args.matrix_height is not None:
        data_selection = args.matrix_height
    elif args.fields is not None:
        data_selection = args.fields

    colorbar_symmetry = 0 if args.symmetry else None
    p = kplt.GenomicVectorArrayPlot(array, keys=data_selection, y_scale=norm, colormap=args.colormap,
                                    colorbar_symmetry=colorbar_symmetry, vmin=args.vmin, vmax=args.vmax,
                                    show_colorbar=args.show_colorbar, replacement_color=args.replacement_color)
    return p, args


def region_parser():
    parser = subplot_parser()
    parser.description = '''Region plot.'''

    parser.add_argument(
        'bed',
        help='''BED or other genomic coordinate file .'''
    )

    parser.add_argument(
        '-f', '--features', dest='features',
        nargs='+',
        help='''(Only for GTF) Plot only the specified feature types (3rd GTF column).'''
    )

    parser.add_argument(
        '-cn', '--color-neutral', dest='color_neutral',
        default='grey',
        help='''Neutral color (no strand information)'''
    )

    parser.add_argument(
        '-cf', '--color-forward', dest='color_forward',
        default='red',
        help='''Forward color (strand +)'''
    )

    parser.add_argument(
        '-cr', '--color-reverse', dest='color_reverse',
        default='blue',
        help='''Reverse color (strand -)'''
    )

    parser.add_argument(
        '-L', '--no-labels', dest='show_labels',
        action='store_false',
        help='''Do not show element labels.'''
    )
    parser.set_defaults(show_labels=True)
    return parser


def region(parameters):
    parser = region_parser()

    args = parser.parse_args(parameters)

    import kaic.plotting as kplt

    p = kplt.GenomicFeatureScorePlot(os.path.expanduser(args.bed), feature_types=args.features,
                                     color_neutral=args.color_neutral, color_forward=args.color_forward,
                                     color_reverse=args.color_reverse, show_labels=args.show_labels)
    return p, args


def line_parser():
    parser = subplot_parser()
    parser.description = '''Line plot.'''

    parser.add_argument(
        'array',
        help='''Array object, e.g. InsulationIndex, DirectionalityIndex, ... .'''
    )

    parser.add_argument(
        '-f', '--fields', dest='fields',
        nargs='+',
        help='''Only plot these fields, otherwise all are plotted.'''
    )

    parser.add_argument(
        '-r', '--range', dest='range',
        nargs=2,
        type=int,
        help='''Range of y-values to plot (<min> <max> inclusive)'''
    )

    parser.add_argument(
        '-v', '--values', dest='values',
        nargs='+',
        type=int,
        help='''Y-values to plot'''
    )

    parser.add_argument(
        '-y', '--ylim', dest='ylim',
        nargs=2,
        type=float,
        help='''Y-axis limits.'''
    )
    return parser


def line(parameters):
    parser = line_parser()
    args = parser.parse_args(parameters)

    import kaic
    import kaic.plotting as kplt

    array = kaic.load(args.array, mode='r')

    data_selection = None
    if args.range is not None:
        data_selection = []
        for i, y in enumerate(array.y_values):
            if args.range[0] <= y <= args.range[1]:
                data_selection.append(array.data_field_names[i])
    elif args.values is not None:
        data_selection = []
        values_set = set(args.values)
        for i, y in enumerate(array.y_values):
            if int(y) in values_set:
                data_selection.append(array.data_field_names[i])
    elif args.fields is not None:
        data_selection = args.fields

    p = kplt.GenomicRegionsPlot(array, attributes=data_selection, ylim=args.ylim)
    return p, args


def bigwig_parser():
    parser = subplot_parser()
    parser.description = '''BigWig plot.'''

    parser.add_argument(
        'bigwig',
        nargs='+',
        help='''BigWig file(s).'''
    )

    parser.add_argument(
        '-n', '--names', dest='names',
        nargs='+',
        help='''Names for each bigWig (must be same length as number of bigWigs).'''
    )

    parser.add_argument(
        '-b', '--bin', dest='bin',
        type=int,
        help='''Bin bigWig values to genomic region of this fixed size.'''
    )

    parser.add_argument(
        '-y', '--ylim', dest='ylim',
        nargs=2,
        type=float,
        help='''Y-axis limits.'''
    )
    return parser


def bigwig(parameters):
    parser = bigwig_parser()
    args = parser.parse_args(parameters)

    import kaic
    import kaic.plotting as kplt

    bigwigs = []
    for file_name in args.bigwig:
        bigwigs.append(kaic.load(file_name, mode='r'))

    p = kplt.BigWigPlot(bigwigs, names=args.names, bin_size=args.bin, ylim=args.ylim)

    return p, args


def gene_parser():
    parser = subplot_parser()
    parser.description = '''Gene plot.'''

    parser.add_argument(
        'genes',
        help='''Genes file (GTF, BED).'''
    )

    parser.add_argument(
        '-f', '--feature-types', dest='feature_types',
        nargs='+',
        default=('exon',),
        help='''Feature types to be plotted. By default only exons are shown.'''
    )

    parser.add_argument(
        '-cn', '--color-neutral', dest='color_neutral',
        default='gray',
        help='''Color for genes without strand information.'''
    )

    parser.add_argument(
        '-cf', '--color-forward', dest='color_forward',
        default='orangered',
        help='''Color for genes on the '+' strand.'''
    )

    parser.add_argument(
        '-cr', '--color-reverse', dest='color_reverse',
        default='darkturquoise',
        help='''Color for genes on the '-' strand.'''
    )

    parser.add_argument(
        '-cs', '--color-score', dest='color_score',
        action='store_true',
        default=False,
        help='''Assign color based on score value.'''
    )

    parser.add_argument(
        '-v', '--vertical-distance', dest='vdist',
        type=float,
        default=0.2,
        help='''Vertical distance between rows of genes in the plot.'''
    )

    parser.add_argument(
        '-b', '--box-height', dest='box_height',
        type=float,
        default=0.5,
        help='''Height of exon boxes.'''
    )

    parser.add_argument(
        '-s', '--font-size', dest='font_size',
        type=float,
        default=9.,
        help='''Font size for gene labels.'''
    )

    parser.add_argument(
        '-a', '--arrow-size', dest='arrow_size',
        type=float,
        default=8.,
        help='''Size of directionality arrows.'''
    )

    parser.add_argument(
        '-l', '--line-width', dest='line_width',
        type=float,
        default=1.,
        help='''Width of the line along the length of the gene.'''
    )

    parser.add_argument(
        '-g', '--group-by', dest='group_by',
        default='transcript_id',
        help='''Group exons according to this attribute (e.g. for
                            plotting multiple transcripts of the same gene.
                            Default: transcript_id'''
    )

    parser.add_argument(
        '-t', '--text-position', dest='text_position',
        default='alternate',
        help='''Position of gene labels. Can be one of 'top' (above gene), 'bottom' (below gene),
                    and 'alternate' (default, alternating between above and below gene).'''
    )

    parser.add_argument(
        '-m', '--min-gene-size', dest='min_gene_size',
        type=int,
        help='''If set, forces the plotted gene to be at least this wide.'''
    )

    parser.add_argument(
        '-L', '--no-labels', dest='show_labels',
        action='store_false',
        help='''Do not plot labels on genes.'''
    )
    parser.set_defaults(show_labels=True)

    parser.add_argument(
        '-C', '--collapse', dest='collapse',
        action='store_true',
        help='''Collapse all genes onto one row.'''
    )
    parser.set_defaults(collapse=False)
    return parser


def gene(parameters):
    parser = gene_parser()
    args = parser.parse_args(parameters)

    import kaic.plotting as kplt

    genes_file = os.path.expanduser(args.genes)

    p = kplt.GenePlot(genes_file, feature_types=args.feature_types,
                      color_neutral=args.color_neutral, color_forward=args.color_forward,
                      color_reverse=args.color_reverse, color_score=args.color_score,
                      vdist=args.vdist, box_height=args.box_height, font_size=args.font_size,
                      arrow_size=args.arrow_size, line_width=args.line_width,
                      group_by=args.group_by, text_position=args.text_position,
                      show_labels=args.show_labels, collapse=args.collapse)

    return p, args
