#!/bin/bash

# remove dupulicate
bcftools norm -d both --threads=4 output/geno/barley/caigebarley__53355variants__807individuals.vcf -O z -o output/geno/barley/caigebarley__53355variants__807individuals_mkdup.vcf

# imputation
beagle gt=output/geno/barley/caigebarley__53355variants__807individuals.vcf out=output/geno/barley/caigebarley__53355variants__807individuals_imp
