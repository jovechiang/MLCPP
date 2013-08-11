function [score] = sample_score(scoretriplet,sigma)

X = lognrnd(scoretriplet,sigma);
score=mean(X);

