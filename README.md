# NASA Hybrid-Electric-Plane
Based on NASA design challenge requirments, this is the code simulating the flight of a hybrid-electrical passenger plane. The purpose is to outline overall feasibility and give a technical basis for design choices.
There should be 3 functions that take as an input the aircraft state (ex, fuel left/used, battery state, altitude, etc) and output the final state. A main function will define the paramaters and call these three functions in order:

takeOff.m - Takeoff analysis

cruise.m - Crusing flight analysis

landing.m - Landing analysis

planeTestMain.m runs the function planeTestFun.m with parameter changes, so that you can see the effects on an upper estimate of flight time. This is useful for evaluating trade-offs and justifying decisions.
