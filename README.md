# InterMod
Modeling plant-pollinator interactions using a Environmental Niche Modeling Framework

Across their geographic ranges, species may interact any place they coocur, or they may interact in only a subset of the places where they coocur. As a result, the actual presence of an interaction may be the result of environmental contexts. This may often be due to the presence or absence of alternative partners (ex. I interact with a low-quality partner only if high-quality partners are not around), or it may actually be due to environmental conditions (as a pollinator, my preferences for or phenological overlap with a certain plant may change as a function of abiotic conditions).

As I'm thinking about this now, I'm interested in using an ENM approach to model plant a pollinator distributions seperately, and then predict their interactions specifically, comparing the model to the null case of "if they interact once, they interact wherever they coocur." I imagine the divergence between these two is going to be pretty fundamentally linked to the degree of interaction specialism of the partners in questions; if pollinators are quite specialized then they probably interact whenever they coocur with the plant clade(s) they specialize. If the species are more generalist, the actual chance of interaction may be more likely to be contextually dependent on the environment. 

Repository structure:
```{bash}

├── README.md. \
├── .gitignore \
├── Analysis \
│   ├── Premise.md Exploratory Data Analysis
├── Data \
│   ├── Outputs \
├── Manuscript \
├── Protocol \

```