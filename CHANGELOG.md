# Changelog

All notable changes to this project will be documented in this file.

## [1.3.0] - 2025-09-29
### Added 
- New exception `AnnDataNoneInGeneralMetadata` to handle cases where required metadata fields contain None or empty values.

### Changed
- General metadata check in `UploadValidator` now checks if any of general metadata or its ontology term ID exists in the `obs` dataframe. For example, the validator will pass if either `tissue` or `tissue_ontology_term_id` exists and non empty. Otherwise, the validator will raise `AnnDataMissingObsColumns`.


## [1.2.0] - 2025-05-30
### Added
- Strict requirement to have dense matrix with embeddings in obsm. Data Frames and sparse matrices will be ignored.


## [1.1.0] - 2025-03-12
### Added
- CLI interface


## [1.0.0] - 2025-03-11
### Added
- `UploadValidator` class to validate AnnData files
- Gene mapping support for *Homo sapiens* and *Mus musculus*
- Custom error handling with `CapMultiException`
- Basic unit tests for core functionality

### Notes
- This is the first public version of the package.
