# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 10:43:39 2018

@author: rmgu

Utilities for creating SNAP GPT graphs from selected operators and executing them.

"""

import os
import tempfile
import subprocess

from xml.etree.cElementTree import Element, SubElement, tostring

EPSG4326 = "EPSG:4326"


class SnapGraph:

    def __init__(self, graph_id="Graph", version="1.0"):
        self.graph = Element("graph", {'id': graph_id})
        SubElement(self.graph, "version").text = version

    def run(self, gpt_path, parallelism="1", cache_size="2G", graph_filename=None):

        is_temp = False
        if not graph_filename:
            graph_filename = tempfile.mkstemp(suffix="_gpt_graph")[1]
            is_temp = True

        self._indent_XML(self.graph)
        with open(graph_filename, "w+b") as graph_file:
            graph_file.write(tostring(self.graph))
            graph_file.flush()

        command = ['"'+gpt_path+'"', '"'+graph_filename+'"', "-q", str(parallelism),
                   "-c", cache_size, "-x"]
        # See https://github.com/senbox-org/snap-engine/issues/80 for reson for
        # LD_LIBRARY_PATH=.
        if os.name == "posix":
            #command.append("LD_LIBRARY_PATH=.")
            pass
        progress = subprocess.Popen(" ".join(command),
                                    shell=True,
                                    universal_newlines=True,
                                    stdout=subprocess.PIPE,
                                    stdin=open(os.devnull),
                                    stderr=subprocess.STDOUT,).stdout
        for line in iter(progress.readline, ''):
            print(line)

        if is_temp:
            try:
                os.remove(graph_filename)
            except PermissionError:
                pass

    def read_op(self, file, format_name=None):

        node_id, node, parameters = self._set_node_boilerplate("Read", None)
        self._set_parameter(parameters, "file", file)
        self._set_parameter(parameters, "formatName", format_name)

        self.graph.append(node)
        return node_id

    def productset_reader_op(self, file_list):

        node_id, node, parameters = self._set_node_boilerplate("ProductSet-Reader", None)
        self._set_parameter(parameters, "fileList", file_list)

        self.graph.append(node)
        return node_id

    def write_op(self, source_product, file, format_name="GeoTIFF"):

        node_id, node, parameters = self._set_node_boilerplate("Write", source_product)
        self._set_parameter(parameters, "file", file)
        self._set_parameter(parameters, "formatName", format_name)

        self.graph.append(node)
        return node_id

    def biophysical_op(self, source_product, compute_lai=False, compute_fapar=False,
                       compute_fcover=False, compute_cab=False, compute_cw=False):

        node_id, node, parameters = self._set_node_boilerplate("BiophysicalOp", source_product)
        self._set_parameter(parameters, "computeLAI", compute_lai)
        self._set_parameter(parameters, "computeFapar", compute_fapar)
        self._set_parameter(parameters, "computeFcover", compute_fcover)
        self._set_parameter(parameters, "computeCab", compute_cab)
        self._set_parameter(parameters, "computeCw", compute_cw)

        self.graph.append(node)
        return node_id

    def resample_op(self, source_product, target_resolution=None, upsampling="Nearest",
                    downsampling="First", flag_downsampling="First",
                    resample_on_pyramid_level=False, reference_band=None, target_width=None,
                    target_height=None):

        node_id, node, parameters = self._set_node_boilerplate("Resample", source_product)
        self._set_parameter(parameters, "referenceBand", reference_band)
        self._set_parameter(parameters, "targetWidth", target_width)
        self._set_parameter(parameters, "targetHeight", target_height)
        self._set_parameter(parameters, "targetResolution", target_resolution)
        self._set_parameter(parameters, "upsampling", upsampling)
        self._set_parameter(parameters, "downsampling", downsampling)
        self._set_parameter(parameters, "flagDownsampling", flag_downsampling)
        self._set_parameter(parameters, "resampleOnPyramidLevels", resample_on_pyramid_level)

        self.graph.append(node)
        return node_id

    def s2_resampling_op(self, source_product, resolution=60, upsampling="Bilinear",
                         downsampling="Mean", flag_downsampling="First",
                         resample_on_pyramid_level=False):
        node_id, node, parameters = self._set_node_boilerplate("S2Resampling", source_product)
        self._set_parameter(parameters, "resolution", resolution)
        self._set_parameter(parameters, "upsampling", upsampling)
        self._set_parameter(parameters, "downsampling", downsampling)
        self._set_parameter(parameters, "flagDownsampling", flag_downsampling)
        self._set_parameter(parameters, "resampleOnPyramidLevels", resample_on_pyramid_level)

        self.graph.append(node)
        return node_id

    def band_maths_op(self, source_product, name, expression, ntype="Float32", description=None,
                      unit=None, no_data_value="NaN"):

        node_id, node, parameters = self._set_node_boilerplate("BandMaths", source_product)
        target_bands = SubElement(parameters, "targetBands")
        target_band = SubElement(target_bands, "targetBand")
        self._set_parameter(target_band, "name", name)
        self._set_parameter(target_band, "type", ntype)
        self._set_parameter(target_band, "expression", expression)
        self._set_parameter(target_band, "description", description)
        self._set_parameter(target_band, "unit", unit)
        self._set_parameter(target_band, "noDataValue", no_data_value)

        self.graph.append(node)
        return node_id

    def reproject_op(self, source_product, crs, resampling="Nearest", no_data_value="NaN",
                     wkt_file=None, reference_pixel_x=None, reference_pixel_y=None, easting=None,
                     northing=None, orientation=None, pixel_size_x=None, pixel_size_y=None,
                     width=None, height=None, tile_size_x=None, tile_size_y=None,
                     orthorectify=False, elevation_model_name=None, include_tie_point_grids=True,
                     add_delta_bands=False):

        node_id, node, parameters = self._set_node_boilerplate("Reproject", source_product)
        self._set_parameter(parameters, "wktFile", wkt_file)
        self._set_parameter(parameters, "crs", crs)
        self._set_parameter(parameters, "resampling", resampling)
        self._set_parameter(parameters, "referencePixelX", reference_pixel_x)
        self._set_parameter(parameters, "referencePixelY", reference_pixel_y)
        self._set_parameter(parameters, "easting", easting)
        self._set_parameter(parameters, "northing", northing)
        self._set_parameter(parameters, "orientation", orientation)
        self._set_parameter(parameters, "pixelSizeX", pixel_size_x)
        self._set_parameter(parameters, "pixelSizeY", pixel_size_y)
        self._set_parameter(parameters, "width", width)
        self._set_parameter(parameters, "height", height)
        self._set_parameter(parameters, "tileSizeX", tile_size_x)
        self._set_parameter(parameters, "tileSizeY", tile_size_y)
        self._set_parameter(parameters, "orthorectify", orthorectify)
        self._set_parameter(parameters, "elevationModelName", elevation_model_name)
        self._set_parameter(parameters, "noDataValue", no_data_value)
        self._set_parameter(parameters, "includeTiePointGrids", include_tie_point_grids)
        self._set_parameter(parameters, "addDeltaBands", add_delta_bands)

        self.graph.append(node)
        return node_id

    def subset_op(self, source_product, source_bands, region=None, geo_region=None,
                  sub_sampling_x=1, sub_sampling_y=1, full_swath=False, tie_point_grid_names=None,
                  copy_metadata=False):

        node_id, node, parameters = self._set_node_boilerplate("Subset", source_product)
        self._set_parameter(parameters, "sourceBands", source_bands)
        self._set_parameter(parameters, "region", region)
        self._set_parameter(parameters, "geoRegion", geo_region)
        self._set_parameter(parameters, "subSamplingX", sub_sampling_x)
        self._set_parameter(parameters, "subSamplingY", sub_sampling_y)
        self._set_parameter(parameters, "fullSwath", full_swath)
        self._set_parameter(parameters, "tiePointGridNames", tie_point_grid_names)
        self._set_parameter(parameters, "copyMetadata", copy_metadata)

        self.graph.append(node)
        return node_id

    def pdu_stitching_op(self, source_product, source_product_paths, target_dir):
        node_id, node, parameters = self._set_node_boilerplate("PduStitching", source_product)
        self._set_parameter(parameters, "sourceProductPaths", source_product_paths)
        self._set_parameter(parameters, "targetDir", target_dir)

        self.graph.append(node)
        return node_id

    def add_elevation_op(self, source_product, dem_name="SRTM 3Sec",
                         dem_resampling_method="BICUBIC_INTERPOLATION",
                         elevation_band_name="elevation"):
        node_id, node, parameters = self._set_node_boilerplate("AddElevation", source_product)
        self._set_parameter(parameters, "demName", dem_name)
        self._set_parameter(parameters, "demResamplingMethod", dem_resampling_method)
        self._set_parameter(parameters, "elevationBandName", elevation_band_name)

        self.graph.append(node)
        return node_id

    def mosaic_op(self, source_product_list, variable_list, condition_list, combine="OR",
                  crs="EPSG:4326", orthorectify=False, elevation_model_name="ACE30",
                  resampling="Nearest", west_bound=-15.0, north_bound=75.0, east_bound=30.0,
                  south_bound=35.0, pixel_size_x=0.05, pixel_size_y=0.05):
        node_id, node, parameters = self._set_node_boilerplate("Mosaic", source_product_list)
        variables = SubElement(parameters, "variables")
        for var in variable_list:
            variable = SubElement(variables, "variable")
            self._set_parameter(variable, "name", var["name"])
            self._set_parameter(variable, "expression", var["expression"])
        conditions = SubElement(parameters, "conditions")
        for con in condition_list:
            condition = SubElement(conditions, "condition")
            self._set_parameter(condition, "name", con["name"])
            self._set_parameter(condition, "expression", con["expression"])
            self._set_parameter(condition, "output", con["output"])
        self._set_parameter(parameters, "combine", combine)
        self._set_parameter(parameters, "crs", crs)
        self._set_parameter(parameters, "orthorectify", orthorectify)
        self._set_parameter(parameters, "elevationModelName", elevation_model_name)
        self._set_parameter(parameters, "resampling", resampling)
        self._set_parameter(parameters, "westBound", west_bound)
        self._set_parameter(parameters, "northBound", north_bound)
        self._set_parameter(parameters, "eastBound", east_bound)
        self._set_parameter(parameters, "southBound", south_bound)
        self._set_parameter(parameters, "pixelSizeX", pixel_size_x)
        self._set_parameter(parameters, "pixelSizeY", pixel_size_y)

        self.graph.append(node)
        return node_id

    def binning_op(self, source_product_list, aggregator_list, output_file, variable_list=[],
                   region=None, start_date_time=None, period_duration=None,
                   time_filter_method="NONE", min_data_hour=None, spatial_resolution=9.28,
                   super_sampling=1, max_distance_on_earth=-1, mask_expr="True",
                   output_type="Product", output_format="GeoTIFF", output_binned_data=False,
                   output_mapped_product=True, metadata_properties_file="./metadata.properties",
                   metadata_template_dir=".", metadata_aggregator_name="NAME"):
        node_id, node, parameters = self._set_node_boilerplate("Binning", source_product_list)
        aggregators = SubElement(parameters, "aggregators")
        for agg in aggregator_list:
            aggregator = SubElement(aggregators, "aggregator")
            self._set_parameter(aggregator, "type", agg["type"])
            self._set_parameter(aggregator, "varName", agg["varName"])
            self._set_parameter(aggregator, "targetName", agg["targetName"])
        variables = SubElement(parameters, "variables")
        for var in variable_list:
            variable = SubElement(variables, "variable")
            self._set_parameter(variable, "name", var["name"])
            self._set_parameter(variable, "expression", var["expression"])
        self._set_parameter(parameters, "region", region)
        self._set_parameter(parameters, "startDateTime", start_date_time)
        self._set_parameter(parameters, "periodDuration", period_duration)
        self._set_parameter(parameters, "timeFilterMethod", time_filter_method)
        self._set_parameter(parameters, "minDataHour", min_data_hour)
        self._set_parameter(parameters, "numRows", int(round(20037.5336/spatial_resolution/2.)*2.))
        self._set_parameter(parameters, "superSampling", super_sampling)
        self._set_parameter(parameters, "maxDistanceOnEarth", max_distance_on_earth)
        self._set_parameter(parameters, "maskExpr", mask_expr)
        self._set_parameter(parameters, "outputType", output_type)
        self._set_parameter(parameters, "outputFile", output_file)
        self._set_parameter(parameters, "outputFormat", output_format)
        self._set_parameter(parameters, "outputBinnedData", output_binned_data)
        self._set_parameter(parameters, "outputMappedProduct", output_mapped_product)
        self._set_parameter(parameters, "metadataPropertiesFile", metadata_properties_file)
        self._set_parameter(parameters, "metadataTemplateDir", metadata_template_dir)
        self._set_parameter(parameters, "metadataAggregatorName", metadata_aggregator_name)

        self.graph.append(node)
        return node_id

    def _indent_XML(self, elem, level=0):
        i = "\n" + level*"  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for elem in elem:
                self._indent_XML(elem, level+1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i

    def _set_node_id(self, operator_name):
        id_num = 1
        for node in self.graph.iter():
            operator = node.find("operator")
            if operator is not None and operator.text == operator_name:
                id_num = id_num + 1

        return operator_name + "(" + str(id_num) + ")"

    def _set_node_boilerplate(self, operator_name, source_product):
        node_id = self._set_node_id(operator_name)
        node = Element("node", {"id": node_id})
        SubElement(node, "operator").text = operator_name
        sources = SubElement(node, "sources")
        if isinstance(source_product, str):
            SubElement(sources, "sourceProduct", {"refid": source_product})
        elif isinstance(source_product, list):
            sources = SubElement(sources, "sourceProducts")
            sources.text = ",".join(source_product)
        parameters = SubElement(node, "parameters")

        return node_id, node, parameters

    def _set_parameter(self, parameters, parameter_name, parameter_value):
        parameter = SubElement(parameters, parameter_name)
        if parameter_value is not None:
            parameter.text = str(parameter_value)
