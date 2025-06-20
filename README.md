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

## How do I try it out?

Can you find me a letter that comes before alpha? The game is not even a single file yet; each of the main components of the game is being designed and tested independently as separate test files.

'''main.py''' is, of course, the aspirational game file. This file will get the bigger updates rolled into it as they reach relatively stable states. You can run it now to see the rough layout of the game board.

The test files contain most of what I'm actively working on. '''PhysicsEngine.py''' is the file for the revamped, comprehensive molecular dynamics simulator. '''cardtestfile.py''' contains the prototyped card and token game. 

The current dependencies to install to build the test files are:
- ursina
- numpy
- chempy
which can be installed using python's '''pip''' package manager.

## What will it look like?

It's still very very early days, so there aren't any previews, renders, models, or anything else that can help show others the pictures in my head of how the game will turn out. However, as I work on the project, I'll try to start taking and sharing screenshots as soon as it starts to show where it's going.

For now, I can say that there are plans for different biomes based on different environments for microscopic life, a system to unlock cards and for collectibility and rare card art, and for the board to reflect environmental context and changes to the board state. 

If you can lend a hand, let me know, but for now, I'm having fun figuring out this insane project! I'm teaching myself everything I need to know (with some help and advice from anyone who's actually educated and studied in the relevant fields), so you won't need a degree to make sense of it.

## Okay, but... how will it actually work? Like what does gameplay look and feel like?

The basic play of a game goes like this.

You select a deck and a level, and enter matchmaking or host/join a game. You are loaded in to a game board in which you are positioned within a translucent bubble- your cell.

You will see icons on the right side of your UI featuring a pink 'nucleotide' icon and a blue 'ATP' icon. You start the game with a supply of nucleotides.

You then select the floating orb in your cell- the nucleus- to open your deck menu and select cards you'd like to use in the early game. You are charged nucleotide points for each card developed this way.

As you exit the menu the nucleus shows the progress of developing cards, and as they are completed they fly to the bottom of your screen to your hand. On hovering the mouse over these cards in your hand, they move up and fan out to allow you to read the full text of each card.

The cards feature a name, unique art, a maximum temperature value and pH limits, and text explaining any effects, costs, and products of the cards. Some cards will be directly playable for effects like a cone of light from the mouse, a beam of UV rays, or a drop-bundle of molecules; whereas others can be played to the deck to work as enzymatic 'factories', soaking up and spitting out the appropriate substrates and products. Some will appear darkened with an ATP cost, which you can choose to pay as often and as much as you like from your ATP storage on your sidebar.

You will play a couple of cards to the board and spawn some early-game simple molecules, and usually focus on building an early economy. This will ensure you a supply of nucleotides for more (and often more point-expensive) cards, ATP to fuel card effects, and the appropriate chemicals on which your strategy or economy relies. When you produce DNA/RNA nucleotides or ATP (by card or by chemical reaction), the molecule pops out into a little glowing pink or cyan (appropriately) sprite, which can be collected manually by passing the mouse over them, or will tick automatically into the sidebar if not interacted with for a short time.

You may then gather up some small molecules into a bundle using the collector tool, shrinking them away to a numbered icon on your sidebar below nucleotides and ATP. Unlike those special 'points' molecules, collected molecules can only be dropped as a bundle- you can't drop them here and there. By default, you can also only drop molecules you've collected up to within the boundaries of your own cell.

Around the mid-game, both players will have started to influence the area between both of you in the centre of the board- gathering chemicals, releasing others, modifying the environment and responding to the situation. You and your opponent will have started to develop cards that define your deck by now, or will be actively engaging an early aggressive tactic. Depending on the map, random events may have now changed the boardstate, providing symmetric effects for both players to adapt to.

Perhaps vacuoles have been dispatched with cargo of threatening chemicals towards the opponent, powered by a membrane-bound Flaggelum card. Or maybe a player with reinforced cell membranes is acidifying the field, damaging the opponents' membrane that slows the movement of the acids towards their delicate economy. Maybe otherwise a flush of rich hydrocarbons led a player to pivot to a strategy that employs otherwise expensive lipids that can be produced from the long chains. Real biochemistry is the inspiration, but the method is very much a tabletop battle of small molecules.

By the end of the game, a competent player will have their win condition in sight and be actively expending their resources per their plan to counter the opponents' defenses and destabilise their nucleotide economy. 

As soon as either player is unable to draw cards, the game is over. The exact way this limit works will be the product of the initial testing and design phase of the prototyped game, simultaneously with detailed balance planning and card/deck designing, in which the state of the game when it reaches players will be worked through. First, I need to feel how the fundamentals of the game feel! Then it'll be easier to compliment the fun and challenge, and not ramp up frustration factors.

## How does it work? Why is it such a mad project?

I'm effectively blending simulation software that's used for research purposes into the dynamics of molecules with a game engine. And I'm doing it in *python*, of all things. I'm considering moving it to another engine/language (perhaps unity, cryengine, or something else based on C++,C#, or Java, since I used those many moons ago), but for now I'm trying to get a working prototype to get the feel for the gameplay and to design the really fun stuff- cards, decks, playstyles, metagame, and art.

I'm currently working on a model of spring-potential bonds, Lennard-Jones non-bonding interactions, and EM forces calculated with Coulomb's law for the simulation. I want this thing to feel as real as possible, with the molecules giving that distinctive haptic feel, the bumps and wobbles of dipole interactions, the sparking migration of electrons, the flow of balance between stochiometric companions, the exciting bubbling of quantity into quality as phases change. If this is all foreign to you, that's okay! Most games are when you start out. Once you learn who the main characters are, what the moves you can make are, what boardstates can look like, all the jargon just becomes another set of names for gameplay elements. I don't think this one will be strictly for the chemistry nerds- though it is definitely one for the nerds. 

I'm not sure if I'll want to sophisticate/modify my physics engine by integrating new/other potentials, or changing the way that they work. There's the option of truncating and splining Lennard-Jones, using other potentials that might better model different states of matter, using logic to select from appropriate potential models, integrating spacial partitioning, and so much more. Let me just say, though, that no wavefunctions or quantum effects will be found here; I'm keeping it fairly classical to make this even slightly feasible, and sticking to the simplest and most widely-used tools wherever possible. I'm not at all sure about including any transition metals, but if I do, it'll probably be Iron and one or two others from the first row of d-block. I've heard this should save me some hair...

## What's it made of?
- game engine: __Ursina__
- __numpy__ for smart, fast arrays
- __chempy__ for stochiometry and data
- considering some __OpenGL__ integration for GPU utilisation. Not sure whether to try doing this through Ursina yet, I'll need to research how easy it'll be to use it for the physical chemistry calculations.

## So what's the goal?

Bugfixing and streamlining this thing is gonna be hell. And then, if it all works out, I'll have a lifetime of angry emails from frustrated chemistry nerds whose favourite reaction doesn't work the way it should in my card game. I'd actually be pretty happy with that, honestly:)

I'm really just looking to see if I can make something that satisfies the craving I have to play with acids and ATP and electrons doing their zappy little things. Not like people do in labs, with vials and stoppers and stirrers and so, so much patience, but like a child with a toy, feeling out the way it works. If there will be any educational value to this game, it'll be that it gets some people more familiar with the exciting world of small-molecule chemical interactions, and excited to dig deeper into the world beyond the game. The goal isn't to be educational, though; it's to play!
