#!/bin/bash
cat $1 | cut -f 1,2 -d "*" --complement | sort > $2
