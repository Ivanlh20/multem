clear all;
clc;
y = get_RandomGenerator_CPU(0, 100000, 1);
figure(1);
hist(y, 100)