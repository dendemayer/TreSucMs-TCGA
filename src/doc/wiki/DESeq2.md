# Multi-factor designs  

Experiments with more than one factor influencing the
counts can be analyzed using design formula that include the additional
variables. In fact, DESeq2 can analyze any possible experimental design that
can be expressed with fixed effects terms (multiple factors, designs with
interactions, designs with continuous variables, splines, and so on are all
possible).

By adding variables to the design, one can control for additional variation in
the counts. For example, if the condition samples are balanced across
experimental batches, by including the batch factor to the design, one can
increase the sensitivity for finding differences due to condition. There are
multiple ways to analyze experiments when the additional variables are of
interest and not just controlling factors (see section on [interactions](interactions)).

The data in the pasilla package have a condition of interest (the column
condition), as well as information on the type of sequencing which was
performed (the column type), as we can see below:

```R
colData(dds)
## DataFrame with 7 rows and 3 columns
##            condition        type sizeFactor
##             <factor>    <factor>  <numeric>
## treated1   treated   single-read   1.635501
## treated2   treated   paired-end    0.761216
## treated3   treated   paired-end    0.832660
## untreated1 untreated single-read   1.138338
## untreated2 untreated single-read   1.793541
## untreated3 untreated paired-end    0.649483
## untreated4 untreated paired-end    0.751600
```
We create a copy of the DESeqDataSet, so that we can rerun the analysis using a
multi-factor design.

```R
ddsMF <- dds
```
We change the levels of type so it only contains letters (numbers, underscore
and period are also allowed in design factor levels). Be careful when changing
level names to use the same order as the current levels.

```R
levels(ddsMF$type)
## [1] "paired-end"  "single-read"
levels(ddsMF$type) <- sub("-.*", "", levels(ddsMF$type))
levels(ddsMF$type)
## [1] "paired" "single"
```
We can account for the different types of sequencing, and get a clearer picture
of the differences attributable to the treatment. As condition is the variable
of interest, we put it at the end of the formula. Thus the results function
will by default pull the condition results unless contrast or name arguments
are specified.

Then we can re-run DESeq:

```R
design(ddsMF) <- formula(~ type + condition)
ddsMF <- DESeq(ddsMF)
```
Again, we access the results using the results function.

```R
resMF <- results(ddsMF)
head(resMF)
## log2 fold change (MLE): condition treated vs untreated 
## Wald test p-value: condition treated vs untreated 
## DataFrame with 6 rows and 6 columns
##               baseMean log2FoldChange     lfcSE       stat    pvalue      padj
##              <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric>
## FBgn0000008   95.14429     -0.0405571  0.220040 -0.1843169 0.8537648  0.949444
## FBgn0000014    1.05652     -0.0835022  2.075676 -0.0402289 0.9679106        NA
## FBgn0000017 4352.55357     -0.2560570  0.112230 -2.2815471 0.0225161  0.130353
## FBgn0000018  418.61048     -0.0646152  0.131349 -0.4919341 0.6227659  0.859351
## FBgn0000024    6.40620      0.3089562  0.755886  0.4087340 0.6827349  0.887742
## FBgn0000032  989.72022     -0.0483792  0.120853 -0.4003139 0.6889253  0.890201
```
It is also possible to retrieve the log2 fold changes, p values and adjusted p
values of variables other than the last one in the design. While in this case,
type is not biologically interesting as it indicates differences across
sequencing protocol, for other hypothetical designs, such as ~genotype +
condition + genotype:condition, we may actually be interested in the difference
in baseline expression across genotype, which is not the last variable in the
design.

In any case, the contrast argument of the function results takes a character
vector of length three: the name of the variable, the name of the factor level
for the numerator of the log2 ratio, and the name of the factor level for the
denominator. The contrast argument can also take other forms, as described in
the help page for results and below

```R
resMFType <- results(ddsMF,
contrast=c("type", "single", "paired"))
head(resMFType)
## log2 fold change (MLE): type single vs paired 
## Wald test p-value: type single vs paired 
## DataFrame with 6 rows and 6 columns
##               baseMean log2FoldChange     lfcSE      stat    pvalue      padj
##              <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
## FBgn0000008   95.14429      -0.262373  0.218505 -1.200767 0.2298414  0.536182
## FBgn0000014    1.05652       3.289885  2.052786  1.602644 0.1090133        NA
## FBgn0000017 4352.55357      -0.100020  0.112091 -0.892310 0.3722268  0.683195
## FBgn0000018  418.61048       0.229049  0.130261  1.758388 0.0786815  0.291789
## FBgn0000024    6.40620       0.306051  0.751286  0.407369 0.6837368  0.880472
## FBgn0000032  989.72022       0.237413  0.120286  1.973744 0.0484108  0.217658

```
If the variable is continuous or an interaction term (see section on
interactions) then the results can be extracted using the name argument to
results, where the nam
