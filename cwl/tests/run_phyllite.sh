#!/bin/bash

CWL="cwl/phyllite.cwl"
YAML="cwl/tests/phyllite_config.yaml"

mkdir -p cwl/tests/test_results
RABIX_ARGS="--basedir cwl/tests/test_results"

rabix $RABIX_ARGS $CWL $YAML
