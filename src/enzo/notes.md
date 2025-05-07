1. Grid has Stars which store local star -> True
2. AllStars seems like have all the stars (but why did they need this? can't understand) -> True

Then, should I use AllStars for sending the entire array or processor-processor pair communication?

pair communication will do


~~- localstarlookupmap will be useful (done)~~


- I gotta work on position update for Stars in StarParticleFinalize

- and also Stars should update Grid stars

-  Grid_UpdateStarParticles is done in StarParticleFindAll

 
- multi pair wise communication failed. (too bad)


~~- background acceleration update needed.~~
- each star particle should have a variable of the rank that they belong

- UpdateParticlePositions should be edited.
- Star particles are updated from particle in the grids.


- In RebuildHierarchy, feedback mode of new particles will be adjusted. 
Here, the particle type will become positive again. (in Grid_StarParticleHandler, it's negative when they're born)
- In MirrorToParticle, it seems like particle data is transferred back to Grids, but why?
this rountine can be improved by using StarMap.. why was it written like this in the first place...


- In Grid_DepositMustRefineParticles, allstars is used to flag a cell for refinement, but not sure if they need all the stars from other processors as well.

~~- in the current setup, all stars will be updated only for each processor. Brocast will be done later on only if needed.~~

- background gravity history is only needed when polynomical approximation is applied. in that case, I can include them within star object. For now, I can just use background acceleration within grid.

- Grid has Star list but is not linked to AllStars I think. that's simple stupid.
- Can't really understand why it's written in this extremely unefficient way.


-  /* Set MetaData->NumberOfParticles and prepare TotalStarParticleCountPrevious
           these are to be used in CommunicationUpdateStarParticleCount
           in StarParticleFinalize */ in StarParticleInitialize
           **it seems like it's not good to delete particles right away.**



**- we should output particle information indpenedently from ABYSS in accordance with Enzo output**

- in abyss, pid to index seems like requiring updated during few-body isn't it? or during some re-ordering.



- maybe we have to do merger first when receiving particles.




- what am I gonna do with ramnents?
- ~~ramnents with zero mass should be assigned with different star type~~
  - ~~currently I removed all the zero mass particles~~


- global_variable->time_step vs EnzoTimeStep. there has to be some fix.


- ~~processor and ID matching within ABYSS I think? (done it was simple fix)~~



- I can directly update positions and velocities to the gird using CurrentGrid and GridParticleIndex (star particle position needed for feedback)

- background vector in grid might have room for improvement

- mass should be mutual updated! in case of accretion.

- Deleting (resetting) Background acceleration in grid might need be modified.

- NumberOfSingleParticle happens to be zero 
- EscapeParticle should be properly treated!


**There are three different stars! ( how stupid it is!!!) Grid->Particle (grid bound), Grid->StarParticle (grid bound), Independent processor-wdie StarParticle (AllStars and such)**

## Grid Particle are the basic/ StarParticle only effective for feedback (AllStars) / Grid->StarParticle is just a mediator.
## Update order: 1. Grid Particle -> Grid Star Particle -> Star Particles in StarParticleInitialize(SFA)
## 2. Grid Particle -> AllStars in StarParticleFinalize for feedback
## Our particles based in AllStars (however all connected).
## this is stupidest system I've ever seen.