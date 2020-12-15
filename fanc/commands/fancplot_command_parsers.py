import argparse
import textwrap


def fancplot_parser():
    usage = '''fancplot [<fancplot global parameters>] <region> [<region> ...]
            --plot <plot type> [<plot parameters>] <plot data file(s)> [...]

            Run fancplot --plot <plot type> -h for help on a specific subplot.\n\nPlot types:\n\n'''

    command_descriptions = dict()
    for name, function in globals().items():
        if name.endswith("_parser") and name != 'fancplot_parser':
            parser = function()
            short_name = name[:-7].replace('_', '-')
            command_descriptions[short_name] = parser.description.split(".")[0]

    command_descriptions.pop('type')
    command_descriptions.pop('subplot')

    max_len = max([len(name) for name in command_descriptions.keys()]) + 4

    usage += "-- Matrix --\n"
    for name in ['triangular', 'square', 'split', 'mirror']:
        padding = ' ' * (max_len - len(name))
        usage += "{}{}{}\n".format(name, padding, command_descriptions.pop(name))

    usage += "\n-- Region --\n"
    for name in ['scores', 'line', 'bar', 'layer', 'gene']:
        padding = ' ' * (max_len - len(name))
        usage += "{}{}{}\n".format(name, padding, command_descriptions.pop(name))

    if len(command_descriptions) > 0:
        usage += "\n-- Other --\n"

    for name in command_descriptions.keys():
        padding = ' ' * (max_len - len(name))
        usage += "{}{}{}\n".format(name, padding, command_descriptions.get(name))

    parser = argparse.ArgumentParser(
        description="fancplot plotting tool for fanc",
        usage=textwrap.dedent(usage)
    )

    parser.add_argument(
        'regions',
        nargs='*',
        default=[],
        help='List of region selectors (<chr>:<start>-<end>) or '
             'files with region information (BED, GTF, ...).'
    )

    parser.add_argument(
        '-o', '--output', dest='output',
        help='Suppresses interactive plotting window and redirects plot '
             'to file. Specify path to file when plotting a single region, '
             'and path to a folder for plotting multiple regions.'
    )

    parser.add_argument(
        '-s', '--script', dest='script',
        help='Use a script file to define plot.'
    )

    parser.add_argument(
        '-p', '--plot', dest='plot',
        action='append',
        help='New plot, type will be chosen automatically by '
             'file type, unless "-t" is provided.'
    )

    parser.add_argument(
        '-n', '--name', dest='name',
        default='',
        help='Plot name to be used as prefix when plotting '
             'multiple regions. Is ignored for single region '
             'and interactive plot.'
    )

    parser.add_argument(
        '--width', dest='width',
        type=int,
        default=4,
        help='Width of the figure in inches. Default: 4'
    )

    parser.add_argument(
        '-w', '--window-size', dest='window_size',
        help='Plotting region size in base pairs. If provided, the '
             'actual size of the given region is ignored and instead '
             'a region <chromosome>:<region center - window size/2> - '
             '<region center + window size/2> will be plotted.'
    )

    parser.add_argument(
        '--invert-x', dest='invert_x',
        action='store_true',
        default=False,
        help='Invert x-axis for this plot'
    )

    parser.add_argument(
        '--tick-locations', dest='tick_locations',
        nargs='+',
        help='Manually define the locations of the tick '
             'labels on the genome axis.'
    )

    parser.add_argument(
        '-V', '--version', dest='print_version',
        action='store_true',
        default=False,
        help='Print version information'
    )

    parser.add_argument(
        '--pdf-text-as-font', dest='pdf_text_as_font',
        action='store_true',
        default=False,
        help='When saving a plot to PDF, save text as a font instead of a path. '
             'This will increase the file size, sometimes by a lot, but it makes the '
             'text in plots editable in vector graphics programs such as Inkscape or '
             'Illustrator.'
    )
    return parser


def type_parser():
    parser = argparse.ArgumentParser(
        description="fancplot subplot identifier",
        add_help=False,
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        'type',
        help='Plot type. See fancplot -h for options.'
    )

    parser.add_argument(
        'data',
        nargs='*',
        help='Data to be plotted in subplot.'
    )

    return parser


def subplot_parser(plot_type='<plot type>'):
    parser = argparse.ArgumentParser(
        prog='fancplot <region> -p {}'.format(plot_type),
        description="fancplot subplot identifier",
    )

    parser.add_argument(
        '--aspect-ratio', dest='aspect_ratio',
        type=float,
        help='Aspect ratio of this panel. Default is determined by figure type (usually 1.0).'
    )

    parser.add_argument(
        '--title', dest='title',
        default='',
        help='Title of this plot.'
    )

    parser.add_argument(
        '--fix-chromosome', dest='fix_chromosome',
        action='store_true',
        default=False,
        help='Fix chromosome identifier for this plot '
             '(add or remove "chr" as required). Use this if'
             'there is a mismatch between the nomenclature '
             'used by different datasets in the figure, specifically '
             'if the chromosome prefix for this dataset does not match '
             'the plot region definition.'
    )

    parser.add_argument(
        '--hide-x', dest='hide_x',
        action='store_true',
        default=False,
        help='Hide x-axis for this plot'
    )

    parser.add_argument(
        '--show-minor-ticks', dest='show_minor_ticks',
        action='store_true',
        default=False,
        help='Show minor ticks on genome axis'
    )

    parser.add_argument(
        '--hide-major-ticks', dest='show_major_ticks',
        action='store_false',
        default=True,
        help='Hide major ticks on genome axis.'
    )

    parser.add_argument(
        '--show-tick-legend', dest='show_tick_legend',
        action='store_true',
        default=False,
        help='Show tick legend with distance between ticks '
             'on genome axis'
    )

    return parser


def triangular_parser():
    parser = subplot_parser('triangular')
    parser.description = 'Triangular Hi-C plot.'

    parser.add_argument(
        'hic',
        help='Hi-C object.'
    )

    parser.add_argument(
        '-vmin', '--minimum-value', dest='vmin',
        type=float,
        help='Minimum value assigned the first '
             'color in the colorbar.'
    )

    parser.add_argument(
        '-vmax', '--maximum-value', dest='vmax',
        type=float,
        help='Maximum value assigned the last color '
             'in the colorbar.'
    )

    parser.add_argument(
        '-m', '--maximum-distance', dest='max_dist',
        help='Maximum distance between two points after '
             'which triangle will be truncated.'
    )

    parser.add_argument(
        '-l', '--log', dest='log',
        action='store_true',
        default=False,
        help='Log-scale colorbar.'
    )

    parser.add_argument(
        '-r', '--range-slider', dest='adjust_range',
        action='store_true',
        default=False,
        help='Add vmax/vmin slider to plot'
    )

    parser.add_argument(
        '-u', '--uncorrected', dest='norm',
        action='store_false',
        default=True,
        help='Plot uncorrected Hi-C matrix values.'
    )

    parser.add_argument(
        '-e', '--observed-expected', dest='oe',
        action='store_true',
        default=False,
        help='Log2-O/E transform matrix values.'
    )

    parser.add_argument(
        '-c', '--colormap', dest='colormap',
        help='Matplotlib colormap'
    )

    parser.add_argument(
        '-y', '--ylabel', dest='ylabel',
        default='',
        help='Label for y axis'
    )

    parser.add_argument(
        '-d', '--default_value', dest='default_value',
        type=float,
        help='Which value to use for missing edge '
             'weights (pixels). Default: 0'
    )

    parser.add_argument(
        '-s', '--colorbar-symmetry', dest='colorbar_symmetry',
        type=float,
        default=None,
        help='Make colorbar symmetrical around this value.'
    )

    parser.add_argument(
        '--weight-field', dest='weight_field',
        help='Which value to use for plotting. '
             'Default: weight'
    )

    parser.add_argument(
        '-C', '--no-colorbar', dest='show_colorbar',
        action='store_false',
        default=True,
        help='Do not show colorbar in plot'
    )
    return parser


def square_parser():
    parser = subplot_parser('square')
    parser.description = 'Square Hi-C plot.'

    parser.add_argument(
        'hic',
        help='Hi-C object.'
    )

    parser.add_argument(
        '-vmin', '--minimum-value', dest='vmin',
        type=float,
        help='Minimum value assigned the first '
             'color in the colorbar.'
    )

    parser.add_argument(
        '-vmax', '--maximum-value', dest='vmax',
        type=float,
        help='Maximum value assigned the last '
             'color in the colorbar.'
    )

    parser.add_argument(
        '-u', '--uncorrected', dest='norm',
        action='store_false',
        default=True,
        help='Plot uncorrected Hi-C matrix values.'
    )

    parser.add_argument(
        '-e', '--observed-expected', dest='oe',
        action='store_true',
        default=False,
        help='Log2-O/E transform matrix values. '
             'Automatically sets colormap to bwr and '
             'makes colorbar symmetrical around 0.'
    )

    parser.add_argument(
        '-l', '--log', dest='log',
        action='store_true',
        default=False,
        help='Log-scale colorbar.'
    )

    parser.add_argument(
        '-r', '--range-slider', dest='adjust_range',
        action='store_true',
        default=False,
        help='Add vmax/vmin slider to plot'
    )

    parser.add_argument(
        '-f', '--flip', dest='flip',
        action='store_true',
        default=False,
        help='Flip matrix upside down'
    )

    parser.add_argument(
        '-c', '--colormap', dest='colormap',
        help='Matplotlib colormap'
    )

    parser.add_argument(
        '-s', '--colorbar-symmetry', dest='colorbar_symmetry',
        type=float,
        default=None,
        help='Make colorbar symmetrical around this value.'
    )

    parser.add_argument(
        '-C', '--no-colorbar', dest='show_colorbar',
        action='store_false',
        default=True,
        help='Do not show colorbar in plot'
    )

    parser.add_argument(
        '--weight-field', dest='weight_field',
        help='Which value to use for plotting. '
             'Default: weight'
    )
    return parser


def split_parser():
    parser = subplot_parser('split')
    parser.description = 'Matrix vs matrix plot'

    parser.add_argument(
        'hic_top',
        help='Top Hi-C object.'
    )

    parser.add_argument(
        'hic_bottom',
        help='Bottom Hi-C object.'
    )

    parser.add_argument(
        '-vmin', '--minimum-value', dest='vmin',
        type=float,
        help='Minimum value assigned the first color '
             'in the colorbar.'
    )

    parser.add_argument(
        '-vmax', '--maximum-value', dest='vmax',
        type=float,
        help='Maximum value assigned the last color '
             'in the colorbar.'
    )

    parser.add_argument(
        '-l', '--log', dest='log',
        action='store_true',
        default=False,
        help='Log-scale colorbar'
    )

    parser.add_argument(
        '-r', '--range-slider', dest='adjust_range',
        action='store_true',
        default=False,
        help='Add vmax/vmin slider to plot'
    )

    parser.add_argument(
        '-e', '--observed-expected', dest='oe',
        action='store_true',
        default=False,
        help='Log2-O/E transform matrix values. '
             'Automatically sets colormap to bwr and '
             'makes colorbar symmetrical around 0.'
    )

    parser.add_argument(
        '--colorbar-symmetry', dest='colorbar_symmetry',
        type=float,
        default=None,
        help='Make colorbar symmetrical around this value.'
    )

    parser.add_argument(
        '-c', '--colormap', dest='colormap',
        help='Matplotlib colormap'
    )

    parser.add_argument(
        '-C', '--no-colorbar', dest='show_colorbar',
        action='store_false',
        default=True,
        help='Do not show colorbar in plot'
    )

    parser.add_argument(
        '-s', '--scale-matrices', dest='scaling',
        action='store_true',
        default=False,
        help='Scale matrix values so that they sum up to the same '
             'number of contacts. Since this is potentially very '
             'time-consuming, it is disabled by default, assuming '
             'that matrices are normalised to the same total number '
             'of valid pairs'
    )
    return parser


def mirror_parser():
    parser = subplot_parser()
    parser.description = '"Mirrored" matrix comparison plot'

    parser.add_argument(
        'hic_upper',
        help='Upper Hi-C object.'
    )

    parser.add_argument(
        'hic_lower',
        help='Lower Hi-C object.'
    )

    parser.add_argument(
        '-uvmin', '--minimum-value-upper', dest='vmin_upper',
        type=float,
        help='Minimum value assigned the first color in the colorbar of upper plot.'
    )

    parser.add_argument(
        '-uvmax', '--maximum-value-upper', dest='vmax_upper',
        type=float,
        help='Maximum value assigned the last color in the colorbar of upper plot.'
    )

    parser.add_argument(
        '-lvmin', '--minimum-value-lower', dest='vmin_lower',
        type=float,
        help='Minimum value assigned the first color in '
             'the colorbar of lower plot.'
    )

    parser.add_argument(
        '-lvmax', '--maximum-value-lower', dest='vmax_lower',
        type=float,
        help='Maximum value assigned the last color in '
             'the colorbar of lower plot.'
    )

    parser.add_argument(
        '-d', '--maximum-distance', dest='max_dist',
        help='Maximum distance between two points after '
             'which triangle will be truncated.'
    )

    parser.add_argument(
        '-ul', '--log-upper', dest='log_upper',
        action='store_true',
        default=False,
        help='Log-scale colorbar for upper plot'
    )

    parser.add_argument(
        '-ll', '--log-lower', dest='log_lower',
        action='store_true',
        default=False,
        help='Log-scale colorbar for lower plot'
    )

    parser.add_argument(
        '-uc', '--colormap-upper', dest='colormap_upper',
        help='Matplotlib colormap for upper plot'
    )

    parser.add_argument(
        '-lc', '--colormap-lower', dest='colormap_lower',
        help='Matplotlib colormap for lower plot'
    )

    parser.add_argument(
        '-ue', '--observed-expected-upper', dest='oe_upper',
        action='store_true',
        default=False,
        help='Log2-O/E transform matrix values of upper matrix. '
             'Automatically sets colormap to bwr and '
             'makes colorbar symmetrical around 0.'
    )

    parser.add_argument(
        '-le', '--observed-expected-lower', dest='oe_lower',
        action='store_true',
        default=False,
        help='Log2-O/E transform matrix values of lower matrix. '
             'Automatically sets colormap to bwr and '
             'makes colorbar symmetrical around 0.'
    )

    parser.add_argument(
        '--colorbar-symmetry-upper', dest='colorbar_symmetry_upper',
        type=float,
        default=None,
        help='Make upper colorbar symmetrical around this value.'
    )

    parser.add_argument(
        '--colorbar-symmetry-lower', dest='colorbar_symmetry_lower',
        type=float,
        default=None,
        help='Make lower colorbar symmetrical around this value.'
    )

    parser.add_argument(
        '-C', '--no-colorbars', dest='show_colorbar',
        action='store_false',
        default=True,
        help='Do not show colorbars in plot'
    )
    return parser


def scores_parser():
    parser = subplot_parser()
    parser.description = 'Region scores plot with parameter dependency.'

    parser.add_argument(
        'scores',
        help='Array object, e.g. InsulationScores, DirectionalityIndexes, ... .'
    )

    parser.add_argument(
        '-vmin', '--minimum-value', dest='vmin',
        type=float,
        help='Minimum value assigned the first color in the colorbar.'
    )

    parser.add_argument(
        '-vmax', '--maximum-value', dest='vmax',
        type=float,
        help='Maximum value assigned the last color in the colorbar.'
    )

    parser.add_argument(
        '-l', '--log', dest='log',
        action='store_true',
        default=False,
        help='Log-transform heatmap values'
    )

    parser.add_argument(
        '-r', '--range', dest='range',
        nargs=2,
        type=int,
        help='Range of y-values to plot (<min> <max> inclusive)'
    )

    parser.add_argument(
        '-p', '--parameters', dest='parameters',
        type=str,
        nargs='+',
        help='List of specific window sizes / parameters to plot.'
    )

    parser.add_argument(
        '-c', '--colormap', dest='colormap',
        default='RdBu_r',
        help='Matplotlib colormap (default: RdBu_r)'
    )

    parser.add_argument(
        '-C', '--no-colorbar', dest='show_colorbar',
        action='store_false',
        default=True,
        help='Do not show colorbar in plot'
    )

    parser.add_argument(
        '-S', '--no-symmetry', dest='symmetry',
        action='store_false',
        default=True,
        help='Do not plot colormap symmetrical around 0.'
    )

    parser.add_argument(
        '-g', '--genomic-format', dest='genomic_format',
        action='store_true',
        default=False,
        help='Use abbreviated genomic formatting for y axis '
             'labels (e.g. 1kb instead of 1000).'
    )

    parser.add_argument(
        '-rc', '--replacement-color', dest='replacement_color',
        default='grey',
        help='Color to replace missing values. Default: grey'
    )
    return parser


def line_parser():
    parser = subplot_parser()
    parser.description = '''Line plot.'''

    parser.add_argument(
        'regions',
        nargs='+',
        help='Region-based file, e.g. BED, GFF, BigWig, ...'
    )

    parser.add_argument(
        '-a', '--attribute', dest='attribute',
        default='score',
        help='Attribute to plot. Default: score'
    )

    parser.add_argument(
        '-b', '--bin-size', dest='bin_size',
        help='Bin size. Specify if you want to bin scores into '
             'larger regions (using weighted mean).'
    )

    parser.add_argument(
        '-l', '--labels', dest='labels',
        nargs='+',
        help='Labels for region datasets'
    )

    parser.add_argument(
        '-c', '--colors', dest='colors',
        nargs='+',
        help='Colors for region datasets'
    )

    parser.add_argument(
        '-f', '--fill', dest='fill',
        action='store_true',
        default=False,
        help='Fill region between the line and the x axis.'
    )

    parser.add_argument(
        '-s', '--line-style', dest='line_style',
        default='mid',
        help='Style of line. Default is "mid": Connect regions only at their'
             'midpoint. Alternatively use "step", where the whole region '
             'is assigned a value, leading to rectangular appearance. '
    )

    parser.add_argument(
        '-y', '--ylim', dest='ylim',
        nargs=2,
        type=float,
        help='''Y-axis limits.'''
    )

    parser.add_argument(
        '--alpha', dest='alpha',
        default=0.5,
        type=float,
        help='Transparency of bars. Value between 0 and 1, default: 0.5'
    )

    parser.add_argument(
        '--legend-location', dest='legend_location',
        default='best',
        help='Location of the legend when providing data labels.'
    )
    return parser


def bar_parser():
    parser = subplot_parser()
    parser.description = '''Bar plot for region scores.'''

    parser.add_argument(
        'regions',
        nargs='+',
        help='Region-based file, e.g. BED, GFF, BigWig, ...'
    )

    parser.add_argument(
        '-a', '--attribute', dest='attribute',
        default='score',
        help='Attribute to plot. Default: score'
    )

    parser.add_argument(
        '-b', '--bin-size', dest='bin_size',
        help='Bin size. Specify if you want to bin scores into '
             'larger regions (using weighted mean).'
    )

    parser.add_argument(
        '--alpha', dest='alpha',
        default=0.5,
        type=float,
        help='Transparency of bars. Value between 0 and 1, default: 0.5'
    )

    parser.add_argument(
        '-l', '--labels', dest='labels',
        nargs='+',
        help='Labels for region datasets'
    )

    parser.add_argument(
        '-y', '--ylim', dest='ylim',
        nargs=2,
        type=float,
        help='''Y-axis limits.'''
    )

    parser.add_argument(
        '-c', '--colors', dest='colors',
        nargs='+',
        help='Colors for region datasets'
    )

    parser.add_argument(
        '--legend-location', dest='legend_location',
        default='best',
        help='Location of the legend when providing data labels.'
    )

    return parser


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
        '-b', '--box-height', dest='box_height',
        type=float,
        default=0.1,
        help='''Height of exon boxes. Default: 0.1'''
    )

    parser.add_argument(
        '-s', '--font-size', dest='font_size',
        type=float,
        default=9.,
        help='''Font size for gene labels. Default: 9'''
    )

    parser.add_argument(
        '-a', '--arrow-size', dest='arrow_size',
        type=float,
        default=8.,
        help='''Size of directionality arrows. Default: 8'''
    )

    parser.add_argument(
        '-l', '--line-width', dest='line_width',
        type=float,
        default=1.,
        help='''Width of the line along the length of the gene. Default: 2.'''
    )

    parser.add_argument(
        '-g', '--group-by', dest='group_by',
        default='transcript_id',
        help='''Group exons according to this attribute (e.g. for
                            plotting multiple transcripts of the same gene.
                            Default: transcript_id'''
    )

    parser.add_argument(
        '--label-field', dest='label_field',
        default='name',
        help='''Which field to use as label.'''
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
        help='''Collapse all genes onto a single row.'''
    )
    parser.set_defaults(collapse=False)

    parser.add_argument(
        '-sq', '--squash', dest='squash',
        action='store_true',
        help='''Merge all exons of each grouping unit together.
                Useful especially when setting --group-by to "gene_id" or "gene_symbol".
                Overlapping genes will still draw on separate rows in contrast to --collapse'''
    )
    parser.set_defaults(squash=False)
    return parser


def layer_parser():
    parser = subplot_parser()
    parser.description = 'Layered feature plot.'

    parser.add_argument(
        'features',
        help='Features file (GFF, BED, Tabix).'
    )

    parser.add_argument(
        '-g', '--grouping', dest='grouping_attribute',
        help='GFF attribute to use for grouping into layers. For BED this is always the name column.'
    )

    parser.add_argument(
        '-i', '--include', dest='include',
        nargs='+',
        help='Include only these groups.'
    )

    parser.add_argument(
        '-e', '--exclude', dest='exclude',
        nargs='+',
        help='Exclude these groups.'
    )

    parser.add_argument(
        '-S', '--no-shadow', dest='shadow',
        action='store_false',
        default=True,
        help='Plot element size "as is" without surrounding box to make it more visible.'
    )

    parser.add_argument(
        '-C', '--collapse', dest='collapse',
        action='store_true',
        default=False,
        help='Plot all elements in one row.'
    )

    parser.add_argument(
        '-w', '--shadow-width', dest='shadow_width',
        type=float,
        default=0.005,
        help='Vertical distance between rows of genes in the plot.'
    )

    parser.add_argument(
        '--no-colors', dest='colors',
        action='store_false',
        default=True,
        help='Do not add color to plots.'
    )

    return parser