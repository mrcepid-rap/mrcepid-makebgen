# mrcepid-filterbcf Developer Readme

## Testing

The test data in this repository uses the output from mrcepid-filterbcf. 
Please refer to the readme in mrcepid-filterbcf in order to generate the input data.

The only requirement is that the data should be in the `test/test_data` directory.

## Files needed to run the applet

In order for the applet to safely decide how to create BGEN chunks (with chunk-ends being in non-exonic regions), we need to create a dictionary of the genome in 
terms of gene coordinates. You can create this dictionary by running the following command from the root of the repository:

```bash
python scripts/gene_dict.py
```

This will create a file called `final_dict_public.json` in the home directory of the repository. You should upload this file to DNA Nexus, and use the file-dxid as part 
of the applet for the `igene_dict` parameter.

## Running this app with additional computational resources

This app has the following entry points:

* main

When running this app, you can override the instance type to be used by
providing the ``systemRequirements`` field to ```/applet-XXXX/run``` or
```/app-XXXX/run```, as follows:

    {
      systemRequirements: {
        "main": {"instanceType": "mem2_hdd2_x2"}
      },
      [...]
    }

See https://documentation.dnanexus.com/developer/api/running-analyses/io-and-run-specifications#run-specification%22%3ERunSpecification%3C/a%3E
in the API documentation for more information about the
available instance types.
