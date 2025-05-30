# Hi! Welcome to the dev page for Chembattle <3

This is my recent passion project that's occupied me lately. In a whiff of inspiration, a vision came to me of a 1v1 RTS-style card-based game, in which users create, develop, and utilise chemicals against each other. 

I've not worked on a dev project in years, and honestly wasn't planning on it, but if it helps me waste less time and paper drawing out molecular diagrams for no good reason, I'll take it. I hope it ends up fun for more people than just me!

## What is Chembattle, really?

It's an odd one. I think the clearest image I can ask you to conjure is something between Starcraft and Magic: The Gathering. 

In the classic 1v1 mode, two players meet across from each other, a board-like space open between them. Each sits behind the protective wall of their own cell membrane, where they will build their strategic resources and tools.

The resources are gained in two ways. 

The first is by the main way that players interact with the game systems; through __cards__. The player can access a pre-prepared deck of cards from their nucleus, spending nucleotides to develop up any card they like for use. Some cards have immediate effects, like spawning specific chemicals in the player's nuclear zone, or allowing the player to direct a beam of photons across the board for a few moments. Others sit on the board as permanent features in the player's nuclear zone, cytoplasm, or on the membranes. These cards mimic the effect of proteins, enabling the enzymatic transformation of chemicals, breakdown of fuels for energy (as ATP), transport chemicals into the field between players, and so much more.

These card effects allow players to then utilise the other resources; the molecules free-floating in the environment between the two players.

What are they utilising these molecules for? Well, to win the game! How? By breaking down your opponents ability to compete, and eventually shut down their nucleus by damaging its essential systems. 

Chemicals, therefore, are the only real way for either player to affect the other, and so you are both forced to try to turn the environment to your advantage and drive home a strategy that leaves your opponent powerless to prevent your victory.

## What will it look like?

It's still very very early days, so there aren't any previews, renders, models, or anything else that can help show others the pictures in my head of how the game will turn out. However, as I work on the project, I'll try to start taking and sharing screenshots as soon as it starts to show where it's going.

For now, I can say that there are plans for different biomes based on different environments for microscopic life, a system to unlock cards and for collectibility and rare card art, and for the board to reflect environmental context and changes to the board state. 

If you can lend a hand, let me know, but for now, I'm having fun figuring out this insane project!

## How does it work? Why is it such a mad project?

I'm effectively blending simulation software that's used for research purposes into the dynamics of molecules with a game engine. And I'm doing it in *python*, of all things. I'm considering moving it to another engine/language (perhaps unity, cryengine, or something else based on C++,C#, or Java, since I used those many moons ago), but for now I'm trying to get a working prototype to get the feel for the gameplay and to design the really fun stuff- cards, decks, playstyles, metagame, and art.

I'm currently working on a model of spring-potential bonds, Lennard-Jones non-bonding interactions, and EM forces calculated with Coulomb's law for the simulation. I want this thing to feel as real as possible, with the molecules giving that distinctive haptic feel, the bumps and wobbles of dipole interactions, the sparking migration of electrons, the flow of balance between stochiometric companions, the exciting bubbling of quantity into quality as phases change. If this is all foreign to you, that's okay! Most games are when you start out. Once you learn who the main characters are, what the moves you can make are, what boardstates can look like, all the jargon just becomes another set of names for gameplay elements. I don't think this one will be strictly for the chemistry nerds- though it is definitely one for the nerds. 

I'm not sure if I'll want to sophisticate/modify my physics engine by integrating new/other potentials, or changing the way that they work. There's the option of truncating and splining Lennard-Jones, using other potentials that might better model different states of matter, using logic to select from appropriate potential models, integrating spacial partitioning, and so much more. 

Bugfixing and streamlining this thing is gonna be hell. And then, if it all works out, I'll have a lifetime of angry emails from frustrated chemistry nerds whose favourite reaction doesn't work the way it should in my card game. I'd actually be pretty happy with that, honestly:)
