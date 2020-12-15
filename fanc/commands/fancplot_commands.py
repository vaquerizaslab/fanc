import os
import fanc
from fanc.config import config
import fanc.plotting as kplt
import fanc.commands.fancplot_command_parsers as parsers


def triangular(parameters):
    parser = parsers.triangular_parser()
    args = parser.parse_args(parameters)

    if args.colormap is None:
        colormap = config.colormap_hic if not args.oe else 'bwr'
    else:
        colormap = args.colormap

    matrix = fanc.load(os.path.expanduser(args.hic), mode='r')
    from fanc.tools.general import str_to_int

    norm = "lin" if not args.log else "log"
    return kplt.HicPlot(matrix, colormap=colormap,
                        max_dist=str_to_int(args.max_dist) if args.max_dist is not None else None,
                        norm=norm, vmin=args.vmin,
                        vmax=args.vmax, show_colorbar=args.show_colorbar,
                        adjust_range=args.adjust_range, oe=args.oe, log=args.oe,
                        colorbar_symmetry=0 if args.oe and args.colorbar_symmetry is None
                                            else args.colorbar_symmetry,
                        ylabel=args.ylabel, weight_field=args.weight_field,
                        default_value=args.default_value, matrix_norm=args.norm), args


def square(parameters):
    parser = parsers.square_parser()

    args = parser.parse_args(parameters)
    norm = "lin" if not args.log else "log"

    if args.colormap is None:
        colormap = config.colormap_hic if not args.oe else 'bwr'
    else:
        colormap = args.colormap

    matrix = fanc.load(os.path.expanduser(args.hic), mode='r')
    return kplt.HicPlot2D(matrix, colormap=colormap, norm=norm, vmin=args.vmin, vmax=args.vmax,
                          show_colorbar=args.show_colorbar, adjust_range=args.adjust_range,
                          oe=args.oe, log=args.oe,
                          colorbar_symmetry=0 if args.oe and args.colorbar_symmetry is None
                                              else args.colorbar_symmetry,
                          flip=args.flip, matrix_norm=args.norm,
                          weight_field=args.weight_field), args


def split(parameters):
    parser = parsers.split_parser()

    args = parser.parse_args(parameters)

    if args.colormap is None:
        colormap = config.colormap_hic if not args.oe else 'bwr'
    else:
        colormap = args.colormap

    matrix_bottom = fanc.load(os.path.expanduser(args.hic_bottom), mode='r')
    matrix_top = fanc.load(os.path.expanduser(args.hic_top), mode='r')

    norm = "lin" if not args.log else "log"
    sp = kplt.HicComparisonPlot2D(matrix_top, matrix_bottom, colormap=colormap,
                                  norm=norm, vmin=args.vmin,
                                  adjust_range=args.adjust_range,
                                  vmax=args.vmax, scale_matrices=args.scaling,
                                  oe=args.oe, log=args.oe,
                                  colorbar_symmetry=0 if args.oe and args.colorbar_symmetry is None
                                                      else args.colorbar_symmetry,
                                  show_colorbar=args.show_colorbar)
    return sp, args


def mirror(parameters):
    parser = parsers.mirror_parser()

    args = parser.parse_args(parameters)

    if args.colormap_lower is None:
        colormap_lower = config.colormap_hic if not args.oe_lower else 'bwr'
    else:
        colormap_lower = args.colormap_lower

    if args.colormap_upper is None:
        colormap_upper = config.colormap_hic if not args.oe_upper else 'bwr'
    else:
        colormap_upper = args.colormap_upper

    matrix_upper = fanc.load(os.path.expanduser(args.hic_upper), mode='r')
    matrix_lower = fanc.load(os.path.expanduser(args.hic_lower), mode='r')

    from fanc.tools.general import str_to_int

    norm_upper = "lin" if not args.log_upper else "log"
    upper_plot = kplt.HicPlot(matrix_upper, colormap=colormap_upper, max_dist=str_to_int(args.max_dist),
                              norm=norm_upper, vmin=args.vmin_upper, vmax=args.vmax_upper,
                              oe=args.oe_upper, log=args.oe_upper,
                              colorbar_symmetry=0 if args.oe_upper and args.colorbar_symmetry_upper is None
                              else args.colorbar_symmetry_upper,
                              show_colorbar=args.show_colorbar, adjust_range=False)
    norm_lower = "lin" if not args.log_lower else "log"
    lower_plot = kplt.HicPlot(matrix_lower, colormap=colormap_lower, max_dist=str_to_int(args.max_dist),
                              norm=norm_lower, vmin=args.vmin_lower, vmax=args.vmax_lower,
                              oe=args.oe_lower, log=args.oe_lower,
                              colorbar_symmetry=0 if args.oe_lower and args.colorbar_symmetry_lower is None
                              else args.colorbar_symmetry_lower,
                              show_colorbar=args.show_colorbar, adjust_range=False)
    vsp = kplt.VerticalSplitPlot(upper_plot, lower_plot)
    return vsp, args


def scores(parameters):
    parser = parsers.scores_parser()

    args = parser.parse_args(parameters)

    array = fanc.load(os.path.expanduser(args.scores), mode='r')
    norm = "linear" if not args.log else "log"

    if args.range is not None:
        data_selection = []
        for i, y in enumerate(array._parameters):
            if args.range[0] <= y <= args.range[1]:
                data_selection.append(y)
    elif args.parameters is not None:
        data_selection = args.parameters
    else:
        data_selection = array._parameters

    colorbar_symmetry = 0 if args.symmetry else None
    p = kplt.GenomicVectorArrayPlot(array, parameters=data_selection, y_scale=norm, colormap=args.colormap,
                                    colorbar_symmetry=colorbar_symmetry, vmin=args.vmin, vmax=args.vmax,
                                    show_colorbar=args.show_colorbar, replacement_color=args.replacement_color,
                                    genomic_format=args.genomic_format)
    return p, args


def line(parameters):
    parser = parsers.line_parser()
    args = parser.parse_args(parameters)

    regions = [fanc.load(file_name) for file_name in args.regions]

    from fanc.tools.general import str_to_int

    attribute = args.attribute
    bin_size = str_to_int(args.bin_size)
    labels = args.labels
    colors = args.colors
    fill = args.fill
    line_style = args.line_style
    ylim = args.ylim
    alpha = args.alpha
    legend_location = args.legend_location

    if labels is not None and len(labels) != len(regions):
        parser.error("Number of labels ({}) must be the same as number "
                     "of datasets ({})".format(len(labels), len(regions)))

    p = kplt.LinePlot(regions, bin_size=bin_size, fill=fill, attribute=attribute, labels=labels,
                      style=line_style, ylim=ylim, colors=colors, legend_location=legend_location,
                      plot_kwargs={'alpha': alpha})
    return p, args


def bar(parameters):
    parser = parsers.bar_parser()
    args = parser.parse_args(parameters)

    regions = [fanc.load(file_name) for file_name in args.regions]
    from fanc.tools.general import str_to_int

    attribute = args.attribute
    bin_size = str_to_int(args.bin_size)
    labels = args.labels
    ylim = args.ylim
    colors = args.colors
    alpha = args.alpha

    legend_location = args.legend_location

    if labels is not None and len(labels) != len(regions):
        parser.error("Number of labels ({}) must be the same as number "
                     "of datasets ({})".format(len(labels), len(regions)))

    p = kplt.BarPlot(regions, attribute=attribute, labels=labels,
                     ylim=ylim, plot_kwargs={'alpha': alpha}, colors=colors,
                     legend_location=legend_location, bin_size=bin_size)
    return p, args


def gene(parameters):
    parser = parsers.gene_parser()
    args = parser.parse_args(parameters)

    genes_file = os.path.expanduser(args.genes)

    p = kplt.GenePlot(genes_file, feature_types=args.feature_types,
                      color_neutral=args.color_neutral, color_forward=args.color_forward,
                      color_reverse=args.color_reverse, color_score=args.color_score,
                      box_height=args.box_height, font_size=args.font_size,
                      arrow_size=args.arrow_size, line_width=args.line_width,
                      group_by=args.group_by,
                      show_labels=args.show_labels, collapse=args.collapse,
                      squash=args.squash, label_field=args.label_field)

    return p, args


def layer(parameters):
    parser = parsers.layer_parser()
    args = parser.parse_args(parameters)

    features = args.features
    grouping_attribute = args.grouping_attribute
    shadow = args.shadow
    include = args.include
    exclude = args.exclude
    shadow_width = args.shadow_width
    collapse = args.collapse
    if args.colors:
        colors = ((1, 'red'), (-1, 'blue'))
    else:
        colors = None

    p = kplt.FeatureLayerPlot(features, gff_grouping_attribute=grouping_attribute, colors=colors,
                              include=include, exclude=exclude,
                              shadow=shadow, shadow_width=shadow_width, collapse=collapse)
    return p, args


# def highlight_parser():
#     parser = subplot_parser()
#     parser.description = '''Highlight regions plot.'''
#
#     parser.add_argument(
#         'regions',
#         help='''BED file or any other Bedtools compatible format.'''
#     )
#
#     parser.add_argument(
#         '-c', '--color',
#         default="grey",
#         help='''Color for shaded regions.'''
#     )
#
#     parser.add_argument(
#         '-a', '--alpha',
#         default=.5,
#         type=float,
#         help='''Alpha transparency, between 0 and 1.'''
#     )
#
#     return parser
#
#
# def highlight(parameters):
#     parser = highlight_parser()
#     args = parser.parse_args(parameters)
#     plot_kwargs = dict(
#         color=args.color,
#         fill=args.color,
#         alpha=args.alpha,
#     )
#     p = kplt.HighlightAnnotation(args.regions, None, None, plot_kwargs=plot_kwargs)
#     return p, None
