# Free Energy Nearest Neighbor (FENN)
FENN is an R package implementing a Von Neumann Entropy based distance metric learning model. The method identifies informative and non-informative features in the data. It scales and projects the data onto a subspace of maximal separation between classes. The corresponding optimal weights on each dimension are identified by optimizing the Von Neumann entropy. The model can be applied to single class (unsupervised) or multi class data, bridging the gap between PCA, LDA and distance metric learning. 


## Installation
To install the package directly from GitHub, you need the devtools library.
```{R}
library(devtools)
install_github("kouroshz/fenn", local = FALSE)
```
## Usage
```
fit <- fenn(y~., data)
p <- predict(fit, newdata)
```

## Example
```{R}
library(fenn)
Iris <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),Sp = rep(c("s","c","v"), rep(50,3)))
train <- sample(1:150, 75)
table(Iris$Sp[train])
z <- fenn(Sp ~ ., Iris, prior = c(1,1,1)/3, subset = train)
predict(z, Iris[-train, ])$class
print(z)
plot(z)
```
# Citation
**A Free Energy Based Approach For Distance Metric Learning**, Sho Inaba, Carl T. Fakhry, Rahul V. Kulkarni, Kourosh Zarringhalam. *The 25th ACM SIGKDD Conference on Knowledge Discovery and Data Mining ACM*, August 4-8, 2019, Anchorage, AK, USA.