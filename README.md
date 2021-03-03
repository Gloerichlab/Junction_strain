# Junction_strain
Python based scripts that can extract properties such as junction length from a skeletonized image of a cell monolayer (over time)

Author: Jooske Monster (Gloerich lab PhD student)

Log:
first upload: 03-03-2021
-Picture_junctional_strain_V1
-Measure_junctionlength_V1
-example_dataset

Picture_junctional_strain_V1:
-input: folder with folders of segmented images (using SeedWaterSegmenter: https://pypi.org/project/SeedWaterSegmenter/)
-output: creates an image of the changes in junction length over time

Measure_junctionlength_V1:
-input: stack of segmented images (using SeedWaterSegmenter: https://pypi.org/project/SeedWaterSegmenter/)
-output: excel files with properties of junctions and cells for selected cell and its neighbors in two selected frames

