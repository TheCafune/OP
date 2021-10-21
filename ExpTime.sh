#!/bin/sh
for var in Erdos992 Advogato,weighted Reality PagesGovernment WikiElec Dmela HepPh Anybeat PagesCompany AstroPh CondMat Gplus Slashdot BlogCatalog Buzznet GooglePlus MathSciNet TwitterFollows Flickr FourSquare IMDB YoutubeSnap ;do
    julia -O3 TestTime.jl $var 5,50,5000
done
