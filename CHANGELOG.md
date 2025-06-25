# genomic-medicine-sweden/nallorefs: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.3.0 - [2025-06-25]

### `Added`

- [#21](https://github.com/genomic-medicine-sweden/nallorefs/pull/21) - Added the option to supply local SNV databases with `--local_echtvar_databases`

### `Changed`

- [#21](https://github.com/genomic-medicine-sweden/nallorefs/pull/21) - Updated reference files with updated repeat definitions and rank models

### `Removed`

### `Fixed`

- [#21](https://github.com/genomic-medicine-sweden/nallorefs/pull/21) - Fixed missing parameters in schema

### Parameters
| Old parameter | New parameter               |
| ------------- | --------------------------- |
|               | `--local_echtvar_databases` |

> [!NOTE]
> Parameter has been updated if both old and new parameter information is present.
> Parameter has been added if just the new parameter information is present.
> Parameter has been removed if new parameter information isn't present.

## 0.2.0 - [2025-04-15]

### `Added`

- [#17](https://github.com/genomic-medicine-sweden/nallorefs/pull/17) - Added new and updated some of the current repeat definitions, for more information see https://github.com/Clinical-Genomics/reference-files/pull/97

### `Changed`

- [#11](https://github.com/genomic-medicine-sweden/nallorefs/pull/11) - Changed version back to 0.2.0dev
- [#16](https://github.com/genomic-medicine-sweden/nallorefs/pull/16) - Changed to `max_sv_size 999999999` in the VEP examples, since there can be variants longer than the previous 248387328.

### `Fixed`

- [#17](https://github.com/genomic-medicine-sweden/nallorefs/pull/17) - Fixed version mismatch in rank models, for more information see: https://github.com/Clinical-Genomics/reference-files/pull/98

## 0.1.1 - [2025-04-07]

### `Fixed`

- [#12](https://github.com/genomic-medicine-sweden/nallorefs/pull/12) - Fixed gaps in the SV rank model frequencies

## 0.1.0 - [2025-04-02]

Initial release of genomic-medicine-sweden/nallorefs, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [#2](https://github.com/genomic-medicine-sweden/nallorefs/pull/2) - Prepare for release
- [#1](https://github.com/genomic-medicine-sweden/nallorefs/pull/1) - Initial version

### `Fixed`

### `Dependencies`

### `Deprecated`
