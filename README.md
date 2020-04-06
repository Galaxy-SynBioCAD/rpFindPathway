# rpFindPathway

Compare a single or collection of rpSBML to another rpSBML and calculate the similarity score. Uses the MIRIAM annotation to find the species cross-references and the EC numbers to compare the reactions. Also uses the BRSynth annotations to compare species using the InChI key writtings.

## Getting Started

This is a docker galaxy tools, and thus, the docker needs to be built locally where Galaxy is installed. 

## Input

Required information:

* **-input_json_path**: (string) Path to the json input
* **-inputTar_path**: (string) Path to the input Tar
* **-output_path**: (string) Path to the output

## Output

* **output**: (string) Path to the output SBOL file

## Dependencies

* Base Docker Image: [brsynth/rpreader](https://hub.docker.com/r/brsynth/rpreader)

## Installing

To build the image using the Dockerfile, use the following command:

```
docker build -t brsynth/rpfindpathways-standalone:dev .
```

### Running the tests

To run the test, untar the test.tar.xz file and run the following command:

```
python run.py -input test/test_rpGlobalScore.tar -input_format tar -input_sbol test/test.sbol -output test/test_output.tar
```

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

v0.1

## Authors

* **Melchior du Lac**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou
* Joan Hérisson

### How to cite rpFindPathways?
