# nf-core/proteomicslfq: Troubleshooting

<!-- TODO nf-core: Change this documentation if these parameters/errors are not relevant for your workflow -->

## Input files not found

If only no file, only one input file , or only read one and not read two is picked up then something is wrong with your input file declaration

1. The path must be enclosed in quotes (`'` or `"`)
2. The path must have at least one `*` wildcard character.

If the pipeline can't find your files then you will get the following error

```bash
ERROR ~ Cannot find any spectra matching: *.mzml
```

Note that if your sample name is "messy" then you have to be very particular with your glob specification. A file name like `L1-1-D-2h_S1_L002_X1_001.mzml` can be difficult enough for a human to read. Specifying `*{1,2}*.mzml` wont work give you what you want Whilst `*{X1,X2}*.mzml` will.

## Data organization

The pipeline can't take a list of multiple input files - it takes a glob expression. If your input files are scattered in different paths then we recommend that you generate a directory with symlinked files.

## Extra resources and getting help

If you still have an issue with running the pipeline then feel free to contact us.
Have a look at the [pipeline website](https://github.com/nf-core/proteomicslfq) to find out how.

If you have problems that are related to Nextflow and not our pipeline then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).
