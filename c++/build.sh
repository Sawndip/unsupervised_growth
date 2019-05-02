 #!/bin/sh
 echo "Build createNetworkForGrowth"
 cd createNetworkForGrowth
 make
 echo "finished"
 cd ../chainGrowth
 echo "Build chainGrowth"
 make
 echo "finished"
 echo "Build continueChainGrowth"
 cd ../continueChainGrowth
 make
 echo "finished"
 echo "Build chainTest"
 cd ../chainTest
 make
 echo "finished"
 
