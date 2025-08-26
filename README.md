# ABySS - Advanced BaYesian Site and Substitution models

A BEAST2 package implementing nonreversible substitution models and novel site-mixture approaches. A linguaPhylo package
for simulating nonreversible substitution is also included.

## Installation instructions
The package is still in development and does not yet have a release.

Before building ensure [OpenJFX 17](https://gluonhq.com/products/javafx/) and 
[Apache Ant](https://ant.apache.org/manual/install.html) are installed, along with the latest version of BEAST2.

1. Clone repository and build package
```
git clone https://github.com/jsaghafifar/ABySS.git
cd ABySS
ant build
```
2. Unzip package in BEAST2 library
```
mkdir ~/.beast/2.7/ABySS
cp abyss-beast/dist/ABySS.v0.0.1.zip ~/.beast/2.7/ABySS
cd ~/.beast/2.7/ABySS
unzip ABySS.v0.0.1.zip
rm ../beauti.properties
```
ABySS will now be ready to run in BEAST2.