# Members Selection

This code takes a cluster already processed with a membership probability estimation method and selects the probability cut that best isolates the true member stars.

Currently the only supported method  is [pyUPMASK][1], but it could be generalized to any method that outputs a per star probability value.

The code assumes that the following columns exist in the input data:

> _x, _y, Plx, pmRA, pmDE, probs_final, BP-RP, Gmag

(the last two are just for plotting) The code also assumes that the frame is rectangular in the coordinates space.


## Probability cut

Once a membership probability estimation method is applied, we still need to select the probability value that separates the data set into the most likely two classes: cluster members and field stars. We refer to this value as "probability cut" or `P_cut`.

The naive selection of `P_cut=0.5` is seldom appropriate, as the data sets analyzed are often heavily imbalanced. This value hence produces a members class with a large degree of field stars contamination.


## Method

The core of the method can be described as follows:

1. Select a probability cut value `P_cut`
2. Separate stars into members (`probs_final>=P_cut`) and field (`probs_final<P_cut`) classes
3. Reject outliers in the members class, using the `Plx, pmRA, pmDE` space. Stars rejected as outliers become field stars
4. Obtain the field density in the `_x, _y` space using the stars classified as field stars
5. Using the field density, estimate the number of members in the cluster area defined by the members
6. Store the (absolute valued) difference between the number of stars in the members class, and the number of members estimated in step 5
7. Repeat steps 1-6 changing the `P_cut` by a small step in a given range
8. Select the `P_cut` that minimizes the difference stored in step 6



[1]: https://github.com/msolpera/pyUPMASK