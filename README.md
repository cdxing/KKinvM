# KKinvM
Invariant mass of phi meson (K+K- pairs) in normal events and mixed events
- `RunAnalyzer.C` is used to load library for `PicoDstAnalyzer2.C` .
- `PicoDstAnalyzer2.C` is used to built Kaon TTree.

# 1. Setup/ Sourcing
Login to BNL RACF computer and get the source code:
```
git clone https://github.com/cdxing/KKinvM.git

```

# 2. Environment of Working on RCF
```
stardev
kinit
aklog
cvs co StRoot/StPicoEvent
```

# 3. Compilation and Execution
Inside the /StRoot/StPicoEvent, do

```
make
```
To run the macro
```
root4star -b -q RunAnalyzer.C+
root4star -b -q PicoAnalyzer2.C+'("raw_PicoDst.root","output_name")'
```
where the "2" is the order of flow

# 4.Updating The Package/Setup
```
git pull
```

inside the /StRoot/StPicoEvent
```
make clean
make
```
