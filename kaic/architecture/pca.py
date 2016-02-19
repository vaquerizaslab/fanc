from kaic.data.genomic import MatrixArchitecturalRegionFeature, Edge
import tables as t
import numpy as np


class HicEdgeCollection(MatrixArchitecturalRegionFeature):
    def __init__(self, hics, additional_fields=None, file_name=None, mode='a', tmpdir=None):
        if additional_fields is not None:
            if not isinstance(additional_fields, dict) and issubclass(additional_fields, t.IsDescription):
                # IsDescription subclass case
                additional_fields = additional_fields.columns
        else:
            additional_fields = dict()

        original_fields = additional_fields.copy()

        self.shared_base_field_names = []
        for i, hic in enumerate(hics):
            field_descriptions = hic._edges.coldescrs
            for name, description in field_descriptions.iteritems():

                if name.startswith("_") or name in {'source', 'sink'}:
                    continue

                for j in xrange(len(hics)):
                    new_name = "%s_%d" % (name, j)
                    if new_name in original_fields:
                        raise ValueError("%s already in hic object, please choose different name." % new_name)
                    if new_name in additional_fields:
                        continue

                    if name not in self.shared_base_field_names:
                        self.shared_base_field_names.append(name)

                    description_args = {'dflt': description.dflt,
                                        'pos': len(additional_fields)}
                    if description.__class__ == t.description.StringCol:
                        description_args['itemsize'] = description.itemsize

                    additional_fields[new_name] = description.__class__(**description_args)

        MatrixArchitecturalRegionFeature.__init__(self, file_name=file_name, mode=mode,
                                                  data_fields=additional_fields,
                                                  regions=hics[0].regions, tmpdir=tmpdir)

        self.hics = hics

    def _calculate(self, *args, **kwargs):
        chromosomes = self.hics[0].chromosomes()

        # step 1: combine all hic objects into one
        #         and calculate variance
        for chr_i, chromosome1 in enumerate(chromosomes):
            for chr_j in xrange(chr_i, len(chromosomes)):
                chromosome2 = chromosomes[chr_j]

                edges = dict()
                for i, hic in enumerate(self.hics):
                    for edge in hic.edge_subset(key=(chromosome1, chromosome2), lazy=False):
                        key = (edge.source, edge.sink)
                        if key not in edges:
                            edges[key] = {}
                        for field in self.shared_base_field_names:
                            if field not in edges[key]:
                                edges[key][field] = [None] * len(self.hics)
                            edges[key][field][i] = getattr(edge, field, None)

                for key, d in edges.iteritems():
                    for field, values in d.iteritems():
                        source = key[0]
                        sink = key[1]

                        d = dict()
                        for i, value in enumerate(values):
                            d["%s_%d" % (field, i)] = value

                        self.add_edge(Edge(source=source, sink=sink, **d), flush=False)

        self.flush()


class HicWeightVariance(HicEdgeCollection):
    def __init__(self, hics, file_name=None, mode='a', tmpdir=None):
        additional_fields = {'var': t.Float32Col()}
        HicEdgeCollection.__init__(self, hics, additional_fields=additional_fields, file_name=file_name,
                                   mode=mode, tmpdir=tmpdir)

    def _calculate(self, *args, **kwargs):
        HicEdgeCollection._calculate(self, *args, **kwargs)

        for edge in self.edges(lazy=True):
            weights = np.zeros(len(self.hics))
            for i in xrange(len(self.hics)):
                weights[i] = getattr(edge, 'weight_' + str(i), 0.0)
            edge.var = np.var(weights)
        self.flush()

        self._edges.cols.var.create_csindex()
