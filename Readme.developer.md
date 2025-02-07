# mrcepid-filterbcf Developer Readme

## Testing

Note for development: this applet has now been set to work with local unit tests, to ensure
robust CI/CD.

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
