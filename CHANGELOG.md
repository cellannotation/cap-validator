# Changelog

All notable changes to this project will be documented in this file.

## Unreleased

## [1.5.0] - 2025-11-24

### Added
- New `find_missing_genes()->pd.DataFrame` function was added to UploadValidator to help users identify missing genes in their AnnData files.
- More gene ids were added to Homo Sapiens Gene Map to cover missing genes reported by users from outdated ensembl releases.

## [1.4.1] - 2025-11-17

### Updated

- Gene Map file for Homo Sapiens was updated to support some outdated ENSEMBL ids from the [Issue](https://github.com/cellannotation/cap-validator/issues/9)


## [1.4.0] - 2025-10-14

### Updated
- Gene Map file for Homo Sapiens was updated from [GENCODE Release 49](https://www.gencodegenes.org/human/release_49.html)
- Gene Map file for Mus Musculus was updated from [GENCODE Release M38](https://www.gencodegenes.org/mouse/release_M38.html)


## [1.3.1] - 2025-09-30
### Changed
- Added gene validation for cases when organisms are specified using ontology term IDs instead of names.

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
