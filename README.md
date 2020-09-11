# SNAP-graph
### Introduction
Simple way to create and execute SNAP GPT graphs.

### Example
Mosaic two Sentinel-3 OL_2_LFR images
```
from snap_graph.snap_graph import SnapGraph

input_images = ["/data/S3A_OL_2_LFR____20200901T233502_20200901T233802_20200902T010531_0179_062_201_3600_LN1_O_NR_002.SEN3",
                "/data/S3A_OL_2_LFR____20200901T233202_20200901T233502_20200902T010530_0179_062_201_3420_LN1_O_NR_002.SEN3"]
mosaic_variable_list = [{"name": "OGVI", "expression": "OGVI"}]
mosaic_condition_list = [{"name": "land", "expression": "LQSF.LAND", "output": False}]

graph = SnapGraph()
read_node_ids = []
for image in input_images:
    read_node_ids.append(graph.read_op(image))
mosaic_id = graph.mosaic_op(read_node_ids, mosaic_variable_list, mosaic_condition_list,
                            crs="EPSG:4326", pixel_size_x=0.003, pixel_size_y=0.003,
                            west_bound=135.710339, north_bound=-21.053497999999998,
                            east_bound=156.52340999999998, south_bound=-45.052335)
graph.write_op(mosaic_id, "/data/mosaic.tif", "/data/graph.xml")
graph.run("gpt", "4", "32G")
```