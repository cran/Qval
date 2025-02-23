#' Perform Q-matrix validation methods
#'
#' @description
#' This function uses generalized Q-matrix validation methods to validate the Q-matrix,
#' including commonly used methods such as GDI (de la Torre, & Chiu, 2016; Najera, Sorrel,
#' & Abad, 2019; Najera et al., 2020), Wald (Ma, & de la Torre, 2020), Hull (Najera et al.,
#' 2021), and MLR-B (Tu et al., 2022). It supports different iteration methods (test
#' level or item level; Najera et al., 2020; Najera et al., 2021; Tu et al., 2022) and
#' can apply various attribute search methods (ESA, SSA, PAA; de la Torre, 2008; Terzi, &
#' de la Torre, 2018). More see details.
#'
#' @section The GDI method:
#' The GDI method (de la Torre & Chiu, 2016), as the first Q-matrix validation method
#' applicable to saturated models, serves as an important foundation for various mainstream
#' Q-matrix validation methods.
#'
#' The method calculates the proportion of variance accounted for (\eqn{PVAF}; @seealso \code{\link[Qval]{get.PVAF}})
#' for all possible q-vectors for each item, selects the q-vector with a \eqn{PVAF} just
#' greater than the cut-off point (or Epsilon, EPS) as the correction result, and the variance
#' \eqn{\zeta^2} is the generalized discriminating index (GDI; de la Torre & Chiu, 2016).
#' Therefore, the GDI method is also considered as a generalized extension of the \eqn{delta}
#' method (de la Torre, 2008), which also takes maximizing discrimination as its basic idea.
#' In the GDI method, \eqn{\zeta^2} is defined as the weighted variance of the correct
#' response probabilities across all mastery patterns, that is:
#' \deqn{
#'  \zeta^2 =
#'  \sum_{l=1}^{2^K} \pi_{l} {(P(X_{pi}=1|\mathbf{\alpha}_{l}) - P_{i}^{mean})}^2
#' }
#' where \eqn{\pi_{l}} represents the prior probability of mastery pattern \eqn{l};
#' \eqn{P_{i}^{mean}=\sum_{k=1}^{K}\pi_{l}P(X_{pi}=1|\mathbf{\alpha}_{l})} is the weighted
#' average of the correct response probabilities across all attribute mastery patterns.
#' When the q-vector is correctly specified, the calculated \eqn{\zeta^2} should be maximized,
#' indicating the maximum discrimination of the item. However, in reality, \eqn{\zeta^2}
#' continues to increase when the q-vector is over-specified, and the more attributes that
#' are over-specified, the larger \eqn{\zeta^2} becomes. The q-vector with all attributes set
#' to 1 (i.e., \eqn{\mathbf{q}_{1:K}}) has the largest \eqn{\zeta^2} (de la Torre, 2016).
#' This is because an increase in attributes in the q-vector leads to an increase in item
#' parameters, resulting in greater differences in correct response probabilities across
#' attribute patterns and, consequently, increased variance. However, this increase in
#' variance is spurious. Therefore, de la Torre et al. calculated \eqn{PVAF = \frac{\zeta^2}{\zeta_{1:K}^2}}
#' to describe the degree to which the discrimination of the current q-vector explains
#' the maximum discrimination. They selected an appropriate \eqn{PVAF} cut-off point to achieve
#' a balance between q-vector fit and parsimony. According to previous studies,
#' the \eqn{PVAF} cut-off point is typically set at 0.95 (Ma & de la Torre, 2020; Najera et al., 2021).
#' Najera et al. (2019) proposed using multinomial logistic regression to predict a more appropriate cut-off point for \eqn{PVAF}. 
#' The cut-off point is denoted as \eqn{eps}, and the predicted regression equation is as follows:
#' 
#' \deqn{ 
#' \log \left( \frac{eps}{1-eps} \right) 
#'    = \text{logit}(eps) 
#'    = -0.405 + 2.867 \cdot IQ + 4.840 \times 10^4 \cdot N - 3.316 \times 10^3 \cdot I
#'  }
#'  Where \eqn{IQ} represents the question quality, calculated as the negative difference between the probability of an examinee 
#'  with all attributes answering the question correctly and the probability of an examinee with no attributes answering the question correctly 
#'  (\eqn{IQ = - \left\{ P\left( \mathbf{1} \right) - \left[ 1 - P\left( \mathbf{0} \right) \right] \right\}}), 
#'  and \eqn{N} and \eqn{I} represent the number of examinees and the number of questions, respectively.
#' 
#' 
#' @section The Wald method:
#' The Wald method (Ma & de la Torre, 2020) combines the Wald test with \eqn{PVAF} to correct
#' the Q-matrix at the item level. Its basic logic is as follows: when correcting item \eqn{i},
#' the single attribute that maximizes the \eqn{PVAF} value is added to a vector with all
#' attributes set to \eqn{\mathbf{0}} (i.e., \eqn{\mathbf{q} = (0, 0, \ldots, 0)}) as a starting point.
#' In subsequent iterations, attributes in this vector are continuously added or
#' removed through the Wald test. The correction process ends when the \eqn{PVAF} exceeds the
#' cut-off point or when no further attribute changes occur. The Wald statistic follows an
#' asymptotic \eqn{\chi^{2}} distribution with a degree of freedom of \eqn{2^{K^\ast} - 1}.
#'
#' The calculation method is as follows:
#' \deqn{
#'    Wald = \left[\mathbf{R} \times P_{i}(\mathbf{\alpha})\right]^{'}
#'    (\mathbf{R} \times \mathbf{V}_{i} \times \mathbf{R})^{-1}
#'    \left[\mathbf{R} \times P_{i}(\mathbf{\alpha})\right]
#' }
#' \eqn{\mathbf{R}} represents the restriction matrix; \eqn{P_{i}(\mathbf{\alpha})} denotes
#' the vector of correct response probabilities for item \eqn{i}; \eqn{\mathbf{V}_i} is the
#' variance-covariance matrix of the correct response probabilities for item \eqn{i}, which
#' can be obtained by multiplying the \eqn{\mathbf{M}_i} matrix (de la Torre, 2011) with the
#' variance-covariance matrix of item parameters \eqn{\mathbf{\Sigma}_i}, i.e.,
#' \eqn{\mathbf{V}_i = \mathbf{M}_i \times \mathbf{\Sigma}_i}. The \eqn{\mathbf{\Sigma}_i} can be
#' derived by inverting the information matrix. Using the the empirical cross-product information
#' matrix (de la Torre, 2011) to calculate \eqn{\mathbf{\Sigma}_i}.
#'
#' \eqn{\mathbf{M}_i} is a \eqn{2^{K^\ast} × 2^{K^\ast}} matrix that represents the relationship between
#' the parameters of item \eqn{i} and the attribute mastery patterns. The rows represent different mastery
#' patterns, while the columns represent different item parameters.
#'
#' @section The Hull method:
#' The Hull method (Najera et al., 2021) addresses the issue of the cut-off point in the GDI
#' method and demonstrates good performance in simulation studies. Najera et al. applied the
#' Hull method for determining the number of factors to retain in exploratory factor analysis
#' (Lorenzo-Seva et al., 2011) to the retention of attribute quantities in the q-vector, specifically
#' for Q-matrix validation. The Hull method aligns with the GDI approach in its philosophy
#' of seeking a balance between fit and parsimony. While GDI relies on a preset, arbitrary
#' cut-off point to determine this balance, the Hull method utilizes the most pronounced elbow
#' in the Hull plot to make this judgment. The the most pronounced elbow is determined using
#' the following formula:
#' \deqn{
#'    st = \frac{(f_k - f_{k-1}) / (np_k - np_{k-1})}{(f_{k+1} - f_k) / (np_{k+1} - np_k)}
#' }
#' where \eqn{f_k} represents the fit-index value (can be \eqn{PVAF} @seealso \code{\link[Qval]{get.PVAF}} or
#' \eqn{R2} @seealso \code{\link[Qval]{get.R2}}) when the q-vector contains \eqn{k} attributes,
#' similarly, \eqn{f_{k-1}} and \eqn{f_{k+1}} represent the fit-index value when the q-vector contains \eqn{k-1}
#' and \eqn{k+1} attributes, respectively. \eqn{{np}_k} denotes the number of parameters when the
#' q-vector has \eqn{k} attributes, which is \eqn{2^k} for a saturated model. Likewise, \eqn{{np}_{k-1}}
#' and \eqn{{np}_{k+1}} represent the number of parameters when the q-vector has \eqn{k-1} and
#' \eqn{k+1} attributes, respectively. The Hull method calculates the \eqn{st} index for all possible q-vectors
#' and retains the q-vector with the maximum \eqn{st} index as the corrected result.
#' Najera et al. (2021) removed any concave points from the Hull plot, and when only the first and
#' last points remained in the plot, the saturated q-vector was selected.
#'
#' @section The MLR-B method:
#' The MLR-B method proposed by Tu et al. (2022) differs from the GDI, Wald and Hull method in that
#' it does not employ \eqn{PVAF}. Instead, it directly uses the marginal probabilities of attribute mastery for
#' subjects to perform multivariate logistic regression on their observed scores. This approach assumes
#' all possible q-vectors and conducts \eqn{2^K-1} regression modelings. After proposing regression equations
#' that exclude any insignificant regression coefficients, it selects the q-vector corresponding to
#' the equation with the minimum \eqn{AIC} fit as the validation result. The performance of this method in both the
#' LCDM and GDM models even surpasses that of the Hull method, making it an efficient and reliable
#' approach for Q-matrix correction.
#' 
#' @section The \eqn{\beta} method:
#' The \eqn{\beta} method (Li & Chen, 2024) addresses the Q-matrix validation problem from the 
#' perspective of signal detection theory. Signal detection theory posits that any stimulus is 
#' a signal embedded in noise, where the signal always overlaps with noise. The \eqn{\beta} method 
#' treats the correct q-vector as the signal and other possible q-vectors as noise. The goal is 
#' to identify the signal from the noise, i.e., to correctly identify the q-vector. For a question 
#' \eqn{i} with the q-vector of the \eqn{c}-th type, the \eqn{\beta} index is computed as follows:
#' 
#' \deqn{
#'    \beta_{jc} = \sum_{l=1}^{2^K} \left| \frac{r_{li}}{n_l} P_{ic}(\mathbf{\alpha_l}) - 
#'                 \left(1 - \frac{r_{li}}{n_l}\right) \left[1 - P_{ic}(\mathbf{\alpha_l})\right] \right|
#'               = \sum_{l=1}^{2^K} \left| \frac{r_{li}}{n_l} - \left[1 - P_{ic}(\mathbf{\alpha_l}) \right] \right|
#'  }
#'  
#' In the formula, \eqn{r_{li}} represents the number of examinees in knowledge state \eqn{l} who correctly 
#' answered item \eqn{i}, while \eqn{n_l} is the total number of examinees in knowledge state \eqn{l}. 
#' \eqn{P_{ic}(\mathbf{\alpha_l})} denotes the probability that an examinee in knowledge state \eqn{l} answers 
#' item \eqn{i} correctly when the q-vector for item \eqn{i} is of the \eqn{c}-th type. In fact, 
#' \eqn{\frac{r_{li}}{n_l}} is the observed probability that an examinee in knowledge state \eqn{l} answers 
#' item \eqn{i} correctly, and \eqn{\beta_{jc}} represents the difference between the actual proportion of 
#' correct answers for item \eqn{i} in each knowledge state and the expected probability of answering the 
#' item incorrectly in that state. Therefore, to some extent, \eqn{\beta_{jc}} can be considered as a measure 
#' of discriminability, and the \eqn{\beta} method posits that the correct q-vector maximizes \eqn{\beta_{jc}}, 
#' i.e.:
#' 
#' \deqn{
#'    \mathbf{q}_i = \arg\max_{\mathbf{q}} \left( \beta_{jc} : \mathbf{q} \in \left\{ \mathbf{q}_{ic}, \, c = 1, 2, \dots, 2^{K} - 1 \right\} \right)
#'  }
#'  
#' Therefore, essentially, \eqn{\beta_{jc}} is an index similar to GDI. Both increase as the number of attributes 
#' in the q-vector increases. Unlike the GDI method, the \eqn{\beta} method does not continue to compute 
#' \eqn{\beta_{jc} / \beta_{j[11...1]}} but instead uses the minimum \eqn{AIC} value to determine whether the attributes 
#' in the q-vector are sufficient. In Package Qval, \link[parallel]{parLapply} will be used to accelerate the \eqn{\beta} method.
#' 
#' Please note that the \eqn{\beta} method has different meanings when applying different search algorithms. 
#' For more details, see section 'Search algorithm' below.
#'  
#' @section Iterative procedure:
#' The iterative procedure that one item modification at a time is item level iteration (\code{"item"}) in (Najera
#' et al., 2020, 2021), while the iterative procedure that the entire Q-matrix is modified at each iteration
#' is test level iteration (\code{"test"}) (Najera et al., 2020; Tu et al., 2022).
#' 
#' The steps of the \code{item} level iterative procedure algorithm are as follows:
#' \describe{
#'    \item{Step1}{Fit the \code{CDM} according to the item responses and the provisional Q-matrix (\eqn{\mathbf{Q}^0}).}
#'    \item{Step2}{Validate the provisional Q-matrix and gain a suggested Q-matrix (\eqn{\mathbf{Q}^1}).}
#'    \item{Step3}{for each item, \eqn{PVAF_{0i}} as the \eqn{PVAF} of the provisional q-vector specified in \eqn{\mathbf{Q}^0},
#'                and \eqn{PVAF_{1i}} as the \eqn{PVAF} of the suggested q-vector in \eqn{\mathbf{Q}^1}.}
#'    \item{Step4}{Calculate all items' \eqn{\Delta PVAF_{i}}, defined as \eqn{\Delta PVAF_{i} = |PVAF_{1i} - PVAF_{0i}|}}
#'    \item{Step5}{Define the hit item as the item with the highest \eqn{\Delta PVAF_{i}}.}
#'    \item{Step6}{Update \eqn{\mathbf{Q}^0} by changing the provisional q-vector by the suggested q-vector of the hit item.}
#'    \item{Step7}{Iterate over Steps 1 to 6 until \eqn{\sum_{i=1}^{I} \Delta PVAF_{i} = 0}}
#' }
#' When the Q-matrix validation method is \code{'MLR-B'}, \code{'Hull'} when \code{criter = 'R2'} or \code{'beta'}, \eqn{PVAF} is not used. 
#' In this case, the criterion for determining which question's index will be replaced is \eqn{AIC}, \eqn{R^2}or \eqn{AIC}, respectively.
#' 
#' \code{test.level = 'test.att'} will use a method called the test-attribute iterative procedure (Najera et al., 2021), which 
#' modifies all items in each iteration while following the principle of minimizing changes in the number of attributes.
#' 
#' The steps of the \code{test} level iterative procedure algorithm are as follows:
#' \describe{
#'    \item{Step1}{Fit the \code{CDM} according to the item responses and the provisional Q-matrix (\eqn{\mathbf{Q}^0}).}
#'    \item{Step2}{Validate the provisional Q-matrix and gain a suggested Q-matrix (\eqn{\mathbf{Q}^1}).}
#'    \item{Step3}{Check whether \eqn{\mathbf{Q}^1 = \mathbf{Q}^0}. If \code{TRUE}, terminate the iterative algorithm.
#'              If \code{FALSE}, Update \eqn{\mathbf{Q}^0} as \eqn{\mathbf{Q}^1}.}
#'    \item{Step4}{Iterate over Steps 1 and 3 until one of conditions as follows is satisfied: 1. \eqn{\mathbf{Q}^1 =
#'                  \mathbf{Q}^0}; 2. Reach the max iteration (\code{maxitr}); 3. \eqn{\mathbf{Q}^1} does not satisfy
#'                 the condition that an attribute is measured by one item at least.}
#' }
#' 
#' @section Search algorithm:
#' Three search algorithms are available: Exhaustive Search Algorithm (ESA), Sequential Search Algorithm (SSA), 
#' and Priority Attribute Algorithm (PAA). 
#' ESA is a brute-force algorithm. When validating the q-vector of a particular item, it traverses all possible 
#' q-vectors and selects the most appropriate one based on the chosen Q-matrix validation method. Since there are 
#' \eqn{2^{K-1}} possible q-vectors with \eqn{K} attributes, ESA requires \eqn{2^{K-1}} searches.
#' 
#' SSA reduces the number of searches by adding one attribute at a time to the q-vector in a stepwise manner. 
#' Therefore, in the worst-case scenario, SSA requires \eqn{K(K-1)/2} searches.
#' The detailed steps are as follows:
#' \describe{
#'    \item{Step 1}{Define an empty q-vector \eqn{\mathbf{q}^0=[00...0]} of length \eqn{K}, 
#'                  where all elements are 0.}
#'    \item{Step 2}{Examine all single-attribute q-vectors, which are those formed by 
#'                  changing one of the 0s in \eqn{\mathbf{q}^0} to 1. 
#'                  According to the criteria of the chosen Q-matrix validation method, 
#'                  select the optimal single-attribute q-vector, denoted as \eqn{\mathbf{q}^1}.}
#'    \item{Step 3}{Examine all two-attribute q-vectors, which are those formed by changing 
#'                  one of the 0s in \eqn{\mathbf{q}^1} to 1. According to the criteria of the 
#'                  chosen Q-matrix validation method, select the optimal two-attribute q-vector, 
#'                  denoted as \eqn{\mathbf{q}^2}.}
#'    \item{Step 4}{Repeat this process until \eqn{\mathbf{q}^K} is found, or the stopping criterion 
#'                  of the chosen Q-matrix validation method is met.}
#' }
#' 
#' PAA is a highly efficient and concise algorithm that evaluates whether each attribute needs to be included in the 
#' q-vector based on the priority of the attributes. @seealso \code{\link[Qval]{get.priority}}. Therefore, even in 
#' the worst-case scenario, PAA only requires \eqn{K} searches. The detailed process is as follows:
#' \describe{
#'    \item{Step 1}{Using the applicable CDM (e.g. the G-DINA model) to estimate the model parameters 
#'                  and obtain the marginal attribute mastery probabilities matrix \eqn{\mathbf{\Lambda}}}
#'    \item{Step 2}{Use LASSO regression to calculate the priority of each attribute in the q-vector for item \eqn{i}}
#'    \item{Step 3}{Check whether each attribute is included in the optimal q-vector based on the attribute 
#'                  priorities from high to low seriatim and output the final suggested q-vector according to the 
#'                  criteria of the chosen Q-matrix validation method.}
#' }
#' 
#' It should be noted that the Wald method proposed by Ma & de la Torre (2020) uses a \code{"stepwise"} search approach. 
#' This approach involves incrementally adding or removing 1 from the q-vector and evaluating the significance of 
#' the change using the Wald test: 
#' 1. If removing a 1 results in non-significance (indicating that the 1 is unnecessary), the 1 is removed from the q-vector; 
#'    otherwise, the q-vector remains unchanged. 
#' 2. If adding a 1 results in significance (indicating that the 1 is necessary), the 1 is added to the q-vector; 
#'    otherwise, the q-vector remains unchanged.
#' The process stops when the q-vector no longer changes or when the PVAF reaches the preset cut-off point (i.e., 0.95).
#' Stepwise are unique search approach of the Wald method, and users should be aware of this. Since stepwise is 
#' inefficient and differs significantly from the extremely high efficiency of PAA, Package \code{Qval} also provides \code{PAA} 
#' for q-vector search in the Wald method. When applying the PAA version of the Wald method, the search still 
#' examines whether each attribute is necessary (by checking if the Wald test reaches significance after adding the attribute) 
#' according to attribute priority. The search stops when no further necessary attributes are found or when the 
#' PVAF reaches the preset cut-off point (i.e., 0.95). The "forward" search approach is another search method 
#' available for the Wald method, which is equivalent to \code{'SSA'}. When \code{'Wald'} uses \code{search.method = 'SSA'}, 
#' it means that the Wald method is employing the forward search approach. Its basic process is the same as \code{'stepwise'}, 
#' except that it does not remove elements from the q-vector. Therefore, the "forward" search approach is essentially SSA.
#' 
#' Please note that, since the \eqn{\beta} method essentially selects q-vectors based on \eqn{AIC}, even without using the iterative process, 
#' the \eqn{\beta} method requires multiple parameter estimations to obtain the AIC values for different q-vectors. 
#' Therefore, the \eqn{\beta} method is more time-consuming and computationally intensive compared to the other methods. 
#' Li and Chen (2024) introduced a specialized search approach for the \eqn{\beta} method, which is referred to as the 
#' \eqn{\beta} search (\code{search.method = 'beta'}). The number of searches required is \eqn{2^{K-2} + K + 1}, and 
#' the specific steps are as follows:
#' \describe{
#'    \item{Step 1}{For item \eqn{i}, sequentially examine the \eqn{\beta} values for each single-attribute q-vector, 
#'                  select the largest \eqn{\beta_{most}} and the smallest \eqn{\beta_{least}}, along with the corresponding 
#'                  attributes \eqn{k_{most}} and \eqn{k_{least}}. (K searches)}
#'    \item{Step 2}{Then, add all possible q-vectors (a total of \eqn{2^K - 1}) containing attribute \eqn{k_{most}} and 
#'                  not containing \eqn{k_{least}} to the search space \eqn{\mathbf{S}_i} (a total of \eqn{2^{K-2}})), and unconditionally 
#'                  add the saturated q-vector \eqn{[11\ldots1]} to \eqn{\mathbf{S}_i} to ensure that it is tested.}
#'    \item{Step 3}{Select the q-vector with the minimum AIC from \eqn{\mathbf{S}_i} as the final output of the \eqn{\beta} 
#'                  method. (The remaining \eqn{2^{K-2} + 1} searches)}
#' }
#' The \code{Qval} package also provides three search methods, ESA, SSA, and PAA, for the \eqn{\beta} method. 
#' When the \eqn{\beta} method applies these three search methods, Q-matrix validation can be completed without 
#' calculating any \eqn{\beta} values, as the \eqn{\beta} method essentially uses \code{AIC} for selecting q-vectors. 
#' For example, when applying ESA, the \eqn{\beta} method does not need to perform Step 1 of the \eqn{\beta} search 
#' and only needs to include all possible q-vectors (a total of \eqn{2^K - 1}) in \eqn{\mathbf{S}_i}, then outputs 
#' the corresponding q-vector based on the minimum \eqn{AIC}. When applying SSA or PAA, the \eqn{\beta} method also 
#' does not require any calculation of \eqn{\beta} values. In this case, the \eqn{\beta} method is consistent 
#' with the Q-matrix validation process described by Chen et al. (2013) using relative fit indices. Therefore, when 
#' the \eqn{\beta} method does not use \eqn{\beta} search, it is equivalent to the method of Chen et al. (2013). 
#' To better implement Chen et al. (2013)'s Q-matrix validation method using relative fit indices, the \code{Qval} 
#' package also provides \eqn{BIC}, \eqn{CAIC}, and \eqn{SABIC} as alternatives to validate q-vectors, in addition 
#' to \eqn{AIC}.
#' 
#' @param Y A required \code{N} × \code{I} matrix or data.frame consisting of the responses of \code{N} individuals
#'          to \code{I} items. Missing values need to be coded as \code{NA}.
#' @param Q A required binary \code{I} × \code{K} containing the attributes not required or required, 0 or 1,
#'            to master the items. The \code{i}th row of the matrix is a binary indicator vector indicating which
#'            attributes are not required (coded by 0) and which attributes are required (coded by 1) to
#'            master item \code{i}.
#' @param CDM.obj An object of class \code{CDM.obj}. When it is not NULL, it enables rapid verification
#'                of the Q-matrix without the need for parameter estimation. @seealso \code{\link[Qval]{CDM}}.
#' @param par.method  Type of mtehod to estimate CDMs' parameters; one out of \code{"EM"}, \code{"BM"}. Default = \code{"EM"}
#'                However, \code{"BM"} is only avaible when \code{method = "GDINA"}.
#' @param mono.constraint Logical indicating whether monotonicity constraints should be fulfilled in estimation.
#'                        Default = \code{TRUE}.
#' @param model Type of model to fit; can be \code{"GDINA"}, \code{"LCDM"}, \code{"DINA"}, \code{"DINO"}
#'              , \code{"ACDM"}, \code{"LLM"}, or \code{"rRUM"}. Default = \code{"GDINA"}.
#'              @seealso \code{\link[Qval]{CDM}}.
#' @param method The methods to validata Q-matrix, can be \code{"GDI"}, \code{"Wald"}, \code{"Hull"}, 
#'               \code{"MLR-B"} and \code{"beta"}. The \code{"model"} must be \code{"GDINA"} when 
#'               \code{method = "Wald"}. Please note that the \eqn{\beta} method has different meanings 
#'               when applying different search algorithms. For more details, see section 'Search algorithm' below.
#'               Default = \code{"GDI"}. See details.
#' @param search.method Character string specifying the search method to use during validation.
#'   \describe{
#'     \item{"ESA"}{for exhaustive search algorithm. Can not for the \code{'Wald'} method.}
#'     \item{"SSA"}{for sequential search algorithm (see de la Torre, 2008; Terzi & de la Torre, 2018). 
#'                  This option can be used when the. It will be equal to \code{"forward"} when \code{method = "Wald"}.}
#'     \item{"PAA"}{for priority attribute algorithm. }
#'     \item{"stepwise"}{only for the \code{"Wald" method}}
#'     \item{"beta"}{only for the \code{"beta" method}}
#'   }
#' @param maxitr Number of max iterations. Default = \code{1}.
#' @param iter.level Can be \code{"item"} level, \code{"test.att"} or \code{"test"} level. Default = \code{"test"} and 
#'                   \code{"test.att"} can not for \code{"Wald"} and \code{"MLR-B"}. See details.
#' @param eps Cut-off points of \eqn{PVAF}, will work when the method is \code{"GDI"} or \code{"Wald"}.
#'            Default = \code{0.95}. See details.
#' @param alpha.level alpha level for the wald test. Default = \code{0.05}
#' @param criter The kind of fit-index value. When \code{method = "Hull"}, it can be \code{R^2} for 
#'               \eqn{R_{McFadden}^2} @seealso \code{\link[Qval]{get.R2}} or \code{'PVAF'} for the proportion of 
#'               variance accounted for (\eqn{PVAF}) @seealso \code{\link[Qval]{get.PVAF}}. When 
#'               \code{method = "beta"}, it can be \code{'AIC'}, \code{'BIC'}, \code{'CAIC'} or \code{'SABIC'}.
#'               Default = \code{"PVAF"} for \code{'Hull'} and default = \code{"AIC"} for \code{'beta'}. See details.
#' @param verbose Logical indicating to print iterative information or not. Default is \code{TRUE}
#'
#' @return
#' An object of class \code{validation} is a \code{list} containing the following components:
#' \item{Q.orig}{The original Q-matrix that maybe contains some mis-specifications and need to be validate.}
#' \item{Q.sug}{The Q-matrix that suggested by certain validation method.}
#' \item{process}{A matrix that contains the modification process of each question during each iteration. 
#'        Each row represents an iteration, and each column corresponds to the q-vector index of the respective 
#'        question. The order of the indices is consistent with the row numbering in the matrix generated by 
#'        the @seealso \code{\link[GDINA]{attributepattern}} function in the \code{GDINA} package.}
#' \item{priority}{An \code{I} × \code{K} matrix that contains the priority of every attribute for
#'                 each item. Only when the \code{search.method} is \code{"PAA"}, the value is availble. See details.}
#' \item{Hull.fit}{A \code{list} containing all the information needed to plot the Hull plot, which is 
#'                 available only when \code{method} = \code{"Hull"}.}
#' \item{iter}{The number of iteration.}
#' \item{time.cost}{The time that CPU cost to finish the function.}
#'
#' @author Haijiang Qin <Haijiang133@outlook.com>
#'
#' @references
#'
#' Chen, J., de la Torre, J., & Zhang, Z. (2013). Relative and Absolute Fit Evaluation in Cognitive Diagnosis Modeling. Journal of Educational Measurement, 50(2), 123-140. DOI: 10.1111/j.1745-3984.2012.00185.x 
#' 
#' de la Torre, J., & Chiu, C. Y. (2016). A General Method of Empirical Q-matrix Validation. Psychometrika, 81(2), 253-273. DOI: 10.1007/s11336-015-9467-8.
#'
#' de la Torre, J. (2008). An Empirically Based Method of Q-Matrix Validation for the DINA Model: Development and Applications. Journal of Education Measurement, 45(4), 343-362. DOI: 10.1111/j.1745-3984.2008.00069.x.
#'
#' Li, J., & Chen, P. (2024). A new Q-matrix validation method based on signal detection theory. British Journal of Mathematical and Statistical Psychology, 00, 1–33. DOI: 10.1111/bmsp.12371
#' 
#' Lorenzo-Seva, U., Timmerman, M. E., & Kiers, H. A. (2011). The Hull method for selecting the number of common factors. Multivariate Behavioral Research, 46, 340–364. DOI: 10.1080/00273171.2011.564527.
#'
#' Ma, W., & de la Torre, J. (2020). An empirical Q-matrix validation method for the sequential generalized DINA model. British Journal of Mathematical and Statistical Psychology, 73(1), 142-163. DOI: 10.1111/bmsp.12156.
#'
#' McFadden, D. (1974). Conditional logit analysis of qualitative choice behavior. In P. Zarembka (Ed.), Frontiers in economics (pp. 105–142). New York, NY: Academic Press.
#'
#' Najera, P., Sorrel, M. A., & Abad, F. J. (2019). Reconsidering Cutoff Points in the General Method of Empirical Q-Matrix Validation. Educational and Psychological Measurement, 79(4), 727-753. DOI: 10.1177/0013164418822700.
#'
#' Najera, P., Sorrel, M. A., de la Torre, J., & Abad, F. J. (2020). Improving Robustness in Q-Matrix Validation Using an Iterative and Dynamic Procedure. Applied Psychological Measurement, 44(6), 431-446. DOI: 10.1177/0146621620909904.
#'
#' Najera, P., Sorrel, M. A., de la Torre, J., & Abad, F. J. (2021). Balancing fit and parsimony to improve Q-matrix validation. British Journal of Mathematical and Statistical Psychology, 74 Suppl 1, 110-130. DOI: 10.1111/bmsp.12228.
#'
#' Terzi, R., & de la Torre, J. (2018). An Iterative Method for Empirically-Based Q-Matrix Validation. International Journal of Assessment Tools in Education, 248-262. DOI: 10.21449/ijate.40719.
#'
#' Tu, D., Chiu, J., Ma, W., Wang, D., Cai, Y., & Ouyang, X. (2022). A multiple logistic regression-based (MLR-B) Q-matrix validation method for cognitive diagnosis models: A confirmatory approach. Behavior Research Methods. DOI: 10.3758/s13428-022-01880-x.
#'
#' @examples
#' ################################################################
#' #                           Example 1                          #
#' #             The GDI method to validate Q-matrix              #
#' ################################################################
#' set.seed(123)
#'
#' library(Qval)
#'
#' ## generate Q-matrix and data
#' K <- 4
#' I <- 20
#' example.Q <- sim.Q(K, I)
#' IQ <- list(
#'   P0 = runif(I, 0.0, 0.2),
#'   P1 = runif(I, 0.8, 1.0)
#' )
#' example.data <- sim.data(Q = example.Q, N = 500, IQ = IQ,
#'                          model = "GDINA", distribute = "horder")
#'
#' ## simulate random mis-specifications
#' example.MQ <- sim.MQ(example.Q, 0.1)
#'
#' \donttest{
#' ## using MMLE/EM to fit CDM model first
#' example.CDM.obj <- CDM(example.data$dat, example.MQ)
#'
#' ## using the fitted CDM.obj to avoid extra parameter estimation.
#' Q.GDI.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "GDI")
#'
#'
#' ## also can validate the Q-matrix directly
#' Q.GDI.obj <- validation(example.data$dat, example.MQ)
#'
#' ## item level iteration
#' Q.GDI.obj <- validation(example.data$dat, example.MQ, method = "GDI",
#'                         iter.level = "item", maxitr = 150)
#'
#' ## search method
#' Q.GDI.obj <- validation(example.data$dat, example.MQ, method = "GDI",
#'                         search.method = "ESA")
#'
#' ## cut-off point
#' Q.GDI.obj <- validation(example.data$dat, example.MQ, method = "GDI",
#'                         eps = 0.90)
#'
#' ## check QRR
#' print(zQRR(example.Q, Q.GDI.obj$Q.sug))
#' }
#'
#'
#'
#' ################################################################
#' #                           Example 2                          #
#' #             The Wald method to validate Q-matrix             #
#' ################################################################
#' set.seed(123)
#'
#' library(Qval)
#'
#' ## generate Q-matrix and data
#' K <- 4
#' I <- 20
#' example.Q <- sim.Q(K, I)
#' IQ <- list(
#'   P0 = runif(I, 0.0, 0.2),
#'   P1 = runif(I, 0.8, 1.0)
#' )
#' example.data <- sim.data(Q = example.Q, N = 500, IQ = IQ, model = "GDINA",
#'                          distribute = "horder")
#'
#' ## simulate random mis-specifications
#' example.MQ <- sim.MQ(example.Q, 0.1)
#'
#' \donttest{
#' ## using MMLE/EM to fit CDM first
#' example.CDM.obj <- CDM(example.data$dat, example.MQ)
#'
#' ## using the fitted CDM.obj to avoid extra parameter estimation.
#' Q.Wald.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "Wald")
#'
#'
#' ## also can validate the Q-matrix directly
#' Q.Wald.obj <- validation(example.data$dat, example.MQ, method = "Wald")
#'
#' ## check QRR
#' print(zQRR(example.Q, Q.Wald.obj$Q.sug))
#' }
#'
#'
#'
#' ################################################################
#' #                           Example 3                          #
#' #             The Hull method to validate Q-matrix             #
#' ################################################################
#' set.seed(123)
#'
#' library(Qval)
#'
#' ## generate Q-matrix and data
#' K <- 4
#' I <- 20
#' example.Q <- sim.Q(K, I)
#' IQ <- list(
#'   P0 = runif(I, 0.0, 0.2),
#'   P1 = runif(I, 0.8, 1.0)
#' )
#' example.data <- sim.data(Q = example.Q, N = 500, IQ = IQ, model = "GDINA",
#'                          distribute = "horder")
#'
#' ## simulate random mis-specifications
#' example.MQ <- sim.MQ(example.Q, 0.1)
#'
#' \donttest{
#' ## using MMLE/EM to fit CDM first
#' example.CDM.obj <- CDM(example.data$dat, example.MQ)
#'
#' ## using the fitted CDM.obj to avoid extra parameter estimation.
#' Q.Hull.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "Hull")
#'
#'
#' ## also can validate the Q-matrix directly
#' Q.Hull.obj <- validation(example.data$dat, example.MQ, method = "Hull")
#'
#' ## change PVAF to R2 as fit-index
#' Q.Hull.obj <- validation(example.data$dat, example.MQ, method = "Hull", criter = "R2")
#'
#' ## check QRR
#' print(zQRR(example.Q, Q.Hull.obj$Q.sug))
#' }
#'
#'
#'
#' ################################################################
#' #                           Example 4                          #
#' #             The MLR-B method to validate Q-matrix            #
#' ################################################################
#' set.seed(123)
#'
#' library(Qval)
#'
#' ## generate Q-matrix and data
#' K <- 4
#' I <- 20
#' example.Q <- sim.Q(K, I)
#' IQ <- list(
#'   P0 = runif(I, 0.0, 0.2),
#'   P1 = runif(I, 0.8, 1.0)
#' )
#' example.data <- sim.data(Q = example.Q, N = 500, IQ = IQ, model = "GDINA",
#'                          distribute = "horder")
#'
#' ## simulate random mis-specifications
#' example.MQ <- sim.MQ(example.Q, 0.1)
#'
#' \donttest{
#' ## using MMLE/EM to fit CDM first
#' example.CDM.obj <- CDM(example.data$dat, example.MQ)
#'
#' ## using the fitted CDM.obj to avoid extra parameter estimation.
#' Q.MLR.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "MLR-B")
#'
#'
#' ## also can validate the Q-matrix directly
#' Q.MLR.obj <- validation(example.data$dat, example.MQ, method  = "MLR-B")
#'
#' ## check QRR
#' print(zQRR(example.Q, Q.MLR.obj$Q.sug))
#' }
#' 
#'
#' @export
#'
validation <- function(Y, Q, 
                       CDM.obj=NULL, par.method="EM", mono.constraint=TRUE, model="GDINA", 
                       method="GDI", search.method="PAA", maxitr=1, iter.level="test",
                       eps=0.95, alpha.level=0.05, criter = NULL, verbose = TRUE){
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)

  if(is.numeric(eps)){
    if(eps > 1 || eps < 0){
      stop("eps only accepts values between 0 and 1 !")
    }
  }else if(is.character(eps)){
    if(eps != "logit"){
      stop("eps can only be logit when it is character !")
    }
  }else{
    stop("eps can only be 'logit' or numeric values between 0 and 1 !")
  }
  
  if(ncol(Q) == 1)
    stop("the number of attributes (K) must be greater than 1 !")

  if(all(method != c("GDI", "Wald", "Hull", "MLR-B", "beta")))
    stop("method must be one of c('GDI', 'Wald', 'Hull', 'MLR-B', 'beta') !")
  
  if(all(search.method != c("ESA", "SSA", "PAA", "stepwise", 'beta')))
    stop("method must be one of c('ESA', 'SSA', 'PAA', 'stepwise', 'beta') !")
  
  if(all(iter.level != c("test", "test.att", "item")))
    stop("method must be one of c('test', 'test.att', 'item') !")
  
  if(iter.level == "test.att" & (method == "Wald" | method == "MLR-B"))
    stop("the iterative of 'test.att' can not for 'Wald' and 'MLR-B' for now !")

  if(search.method == "ESA")
    if(method == "Wald")
      stop("Wald can not work with ", search.method, " !")

  if(search.method == "stepwise")
    if(method != "Wald")
      stop(search.method, " is for Wald method only !")
  
  if(search.method == "beta")
    if(method != "beta")
      stop(search.method, " is for \u03B2 method only !")
  
  if(is.null(criter)){
    if(method == "Hull"){
      criter <- "PVAF"
    }else if(method == "beta"){
      criter <- "AIC"
    }
  }
  
  if(method == "beta")
    if(all(criter != c("AIC", "BIC", "CAIC", "SABIC")))
      stop("criter must be one of c('AIC', 'BIC', 'CAIC', 'SABIC') when method = 'beta' !")
  
  if(method == "Hull")
    if(all(criter != c("PVAF", "R2")))
      stop("criter must be one of c('PVAF', 'R2') when method = 'Hull' !")

  if(method == "Wald" && model != "GDINA")
    stop("model must be GDINA when method are Wald !")

  method.print <- method
  if(method == "beta"){
    method.print <- "\u03B2"
  }
  search.method.print <- search.method
  if(search.method == "beta"){
    search.method.print <- "\u03B2"
  }

  if(verbose){
    cat(method.print, " method with ", search.method.print, " in ", iter.level, " level iteration ...\n")
  }

  time.cost <- system.time({
    if(method == "GDI"){
      Qval.obj <- correctQ.GDI(Y, Q, CDM.obj, method=par.method, mono.constraint=mono.constraint, model=model, 
                               search.method=search.method, maxitr=maxitr, iter.level=iter.level, eps=eps, verbose=verbose)
    }else if(method == "Wald"){
      Qval.obj <- correctQ.Wald(Y=Y, Q=Q, CDM.obj=CDM.obj, mono.constraint=mono.constraint, search.method=search.method, 
                                iter.level=iter.level, maxitr=maxitr, eps=eps, alpha.level=alpha.level, verbose=verbose)
    }else if(method == "Hull"){
      Qval.obj <- correctQ.Hull(Y, Q, CDM.obj, method=par.method, mono.constraint=mono.constraint, model=model, 
                                search.method=search.method, maxitr=maxitr, iter.level=iter.level, criter=criter, verbose=verbose)
    }else if(method == "MLR-B"){
      Qval.obj <- correctQ.MLR(Y, Q, CDM.obj, method=par.method, mono.constraint=mono.constraint, model=model, 
                               search.method=search.method, iter.level=iter.level, maxitr=maxitr, verbose=verbose)
    }else if(method == "beta"){
      Qval.obj <- correctQ.beta(Y, Q, CDM.obj, method=par.method, mono.constraint=mono.constraint, model=model, 
                               search.method=search.method, maxitr=maxitr, iter.level=iter.level, criter=criter, verbose=verbose)
    }
      
  })
  
  res <- list(Q.orig = Qval.obj$Q.original, Q.sug = Qval.obj$Q.sug,
              process = Qval.obj$process, priority = Qval.obj$priority, Hull.fit = Qval.obj$Hull.fit, 
              iter = Qval.obj$iter, time.cost = time.cost[1])
  class(res) <- "validation"

  return(res)
}



