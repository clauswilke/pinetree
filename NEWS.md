# Changelog

## Pinetree 0.5.0

- Full release for tRNA dynamics feature.

## Pinetree 0.4.1

- Adds tRNA dynamics as an experimental/beta feature.
- Fixes bug where RBS propensities are calculated incorrectly when manually adding transcripts to simulations.
- Makes all polymerase speed constants doubles/floats instead of ints.

## Pinetree 0.4.0

- Pinetree can now simulate circular genome expression.
- Updates the catch.hpp dependency to a version that compiles on Apple M1 devices. 

## Pinetree 0.3.0

- Support for site-specific RNase binding constants.
- Overlapping genes and genes that overlap with ribosome binding sites are now fully supported.
- Bug fixes. 

## Pinetree 0.2.0

- It is now possible to simulate transcript degradation from internal and external cleavage sites.
- Fixed transcript abundances for simulating translation without transcription.
- Minor bug fixes and more informative error messages.

## Pinetree 0.1.4

- More informative error message message when installation fails because cmake is not installed.

## Pinetree 0.1.0

First initial release.
