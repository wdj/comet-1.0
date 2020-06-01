#!/bin/bash
#------------------------------------------------------------------------------

./settings_advisor --help

./settings_advisor --num_node 6000 --platform Titan \
                   --num_vector 28342758 --num_field 882 \
                   --metric_type ccc --num_way 2

./settings_advisor --num_node 6000 --platform Titan \
                   --num_vector 28342758 --num_field 882 \
                   --metric_type czekanowski --num_way 2

./settings_advisor --num_node 18688 --platform Titan \
                   --num_vector $(( 1000 * 1000 )) --num_field 82 \
                   --metric_type ccc --num_way 3 \
                   --sparse yes








#------------------------------------------------------------------------------
