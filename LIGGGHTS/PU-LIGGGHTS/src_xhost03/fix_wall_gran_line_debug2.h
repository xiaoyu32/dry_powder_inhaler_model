
            double pointOnSegment[3];
            pointOnSegment[0] = mySegment.center[0] + mySegment.segmentParameter * mySegment.orientation[0];
            pointOnSegment[1] = mySegment.center[1] + mySegment.segmentParameter * mySegment.orientation[1];
            pointOnSegment[2] = mySegment.center[2] + mySegment.segmentParameter * mySegment.orientation[2];

            printf("orientation: (%g %g %g), length: %g cylRadius: %g, x: (%g %g %g), segmentParameter: %g pointOnSegment: (%g %g %g), delta (%g %g %g), deltan: %g \n", 
                    mySegment.orientation[0],mySegment.orientation[1],mySegment.orientation[2],
                    mySegment.length, mySegment.radius,
                    mySegment.center[0],mySegment.center[1],mySegment.center[2],
                    mySegment.segmentParameter,
                    pointOnSegment[0],pointOnSegment[1],pointOnSegment[2],
                    delta[0],delta[1],delta[2],deltan
                  );
