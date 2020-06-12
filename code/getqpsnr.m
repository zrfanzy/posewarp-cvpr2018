listing = dir('resultsimg/*.png');
allscore = 0.0;
for i=1:length(listing)
    im1 = imread('resultsimg/' + listing(i).name);
    im2 = imread('groundtruth/' + listing(i).name);
    allscore = allscore + qpsnr(im1, im2)
end
allscore/length(listing)