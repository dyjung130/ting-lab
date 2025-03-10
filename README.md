# ting-lab

EEG and fMRI codes

General Setup (for MacOS, especially.. :( ).
- Mac OS asked lots of stuff... the biggest problem was the "unverified" .mex.. codes from FieldTrip. You can resolve this by checking the following website: https://www.fieldtriptoolbox.org/faq/matlab/mex_osx/
```
sudo xattr -r -d com.apple.quarantine LOCATION_OF_FIELDTRIP
sudo find LOCATION_OF_FIELDTRIP -name \*.mexmaci64 -exec spctl --add {} \
```
