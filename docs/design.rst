Design and Implementation
=========================

Stochastic Model Algorithm
-------------------------------

Pinetree employs the Gillespie Stochastic Model Algorithm (SSA) \cite{Gillespie1977} to model all molecular interactions involved in gene expression, including the movement of individual polymerases on DNA and ribosomes on mRNA. In the SSA, we consider all molecular interactions as reactions, extending from a parent `Reaction` class. The `Model` class contains the core components of the SSA. The Gillespie SSA is defined as follows:


1. **Initialize**. Begin with time :math:`t = 0`. Specify free species and assign them a copy number. Specify a set of reactions and reaction rate constants involving the free species.
2. **Compute propensities**. For all reactions in the system compute the propensity of the reaction occurring. Sum over all propensities for all reactions to get the total propensity. 
3. **Generate random numbers**. Based on the individual propensities computed in step 2, randomly select a reaction to occur from all reactions in the system. Using the sum of all propensities, :math:`\tau`, compute the time that the next reaction is expected to occur.
4. **Execute**. Execute the reaction from step 4, and advance the time :math:`t` by :math:`\tau`.
5. **Iterate**. Repeat steps 2--5 until :math:`t` exceeds some predefined end time of the simulation.


Pinetree follows these steps, while defining several specialized reactions that provide single-molecule tracking of polymerases and ribosomes along DNA and mRNA. We will now describe these specialized reactions, and each step of the SSA, in more detail.

Pooled species-level reactions
------------------------------

An abstract `Reaction` class is the parent class of all reactions in the SSA. `Model` maintains a vector of all reactions in the SSA, and they must inherit from `Reaction` (Fig.~\ref{fig:bridge}). To comply with the SSA framework, all reactions must be capable of three actions: calculating their propensity (Step 2) and executing (Step 4). The `SpeciesReaction` implements a standard SSA reaction. `SpeciesReaction` defines a set of products, a set of reactants, and a rate constant, and computes its propensity from the copy numbers of interacting species and a rate constant. The SSA requires a mesoscopic rate constant, which differs from a macroscopic rate constant typical of deterministic models in that it depends on the reaction volume. `SpeciesReaction` converts macroscopic rate constants to mesoscopic rate constants internally. Thus the same rate constants from deterministic models can be used directly to parameterize Pinetree simulations. 

Upon execution, `SpeciesReaction` increments its product counts and decrements the reactant counts. Pinetree only supports coefficients of one for reactants. Product coefficients can be any value greater than zero. `SpeciesReaction` objects themselves do not maintain counts of molecular species. All molecular species in Pinetree simulation run are tracked by a single instance of `SpeciesTracker`. `SpeciesTracker` maintains and species-to-`SpeciesReaction` maps. This allows `Model` to cache propensities, and only update the propensities of reactions whose products or reactants have just changed.

Individual molecule-level reactions
-----------------------------------

At the single-molecule level, \texttt{Polymer} objects track individual \texttt{Polymerases} as they move along the polymer, detecting collisions and determining when the \texttt{Polymerase} should leave the \texttt{Polymer}. Tracking RNA polymerases on a genome and ribosomes on mRNA transcripts share enough similarities that most of transcription and translation logic is defined in a generic \texttt{Polymer} class. In practice, a Pinetree simulation run uses the child classes \texttt{Genome} and  \texttt{Transcript}, described in later sections. Ribosomes and RNA polymerases differ only in the definition of their member variables (e.g. footprint size, movement rate, step size, and binding interactions) and are thus represented only by the \texttt{Polymerase} class. Likewise, promoters and ribosome binding sites are represented by a single \texttt{Promoter} class and terminators and stop codons by a single \texttt{Terminator} class. The remainder of this section describes the \texttt{Polymer}, \texttt{Polymerase}, \texttt{Promoter}, and \texttt{Terminator} classes generically.

Pinetree defines two specialized reactions that differ from \texttt{SpeciesReaction} to handle single-molecule tracking. The first are called \texttt{Bind} reactions. \texttt{Bind} reactions coordinate the transfer of a \texttt{Polymerase}s from the pooled species level to an individual \texttt{Polymer} object. \texttt{Bind} reactions are specific for each type \texttt{Promoter}-\texttt{Polymerase} binding interaction. \texttt{Bind} treats all exposed \texttt{Promoters} as a pooled species. Upon execution of the \texttt{Bind} reaction, \texttt{Bind} randomly selects a \texttt{Polymer} that contains open \texttt{Promoter} with which to bind. It constructs a new \texttt{Polymerase}, then instructs the \texttt{Polymer} to bind the \texttt{Polymerase} to an open \texttt{Promoter}. \texttt{SpeciesTracker} maintains \texttt{Promoter}-to-\texttt{Polymer} maps. This map allows all components of Pinetree to quickly look up which \texttt{Polymer}s contain a given \texttt{Promoter}. 

The second type of specialized single-molecule reactions are called \texttt{Bridge} reactions. \texttt{Bridge} reactions are thin wrappers around individual \texttt{Polymer}s. Every individually modeled \texttt{Polymer} has an associated \texttt{Bridge} reaction. The \texttt{Bridge} reaction requests a total propensity value from its \texttt{Polymer}. Upon execution, it instructs the \texttt{Polymer} to move one of its \texttt{Polymerase}s.

`Polymer`
^^^^^^^^^

Each \texttt{Polymer} object maintains a vector of \texttt{Polymerase} objects attached to the \texttt{Polymer}, ordered by position. Each \texttt{Polymerase} defines a fixed speed, in base pairs per second, of movement. The \texttt{Polymer} then scales these movement rates and stores them. Each \texttt{Polymer} maintains a vector of scaling factors representing every base in the \texttt{Polymer}. To compute the scaled movement rates for all \texttt{Polymerase} objects, the movement rate of \texttt{Polymerase} is multiplied by a scaling factor corresponding to the position of the \texttt{Polymer}. These scaled movement rates correspond exactly to SSA propensities. The sum of these propensities represent the overall propensity of the \texttt{Polymer}, i.e., a value proportional to the probability of any single \texttt{Polymerase} moving along the \texttt{Polymer}. This overall propensity is reported to \texttt{Model}, which tracks all \texttt{Polymers} via \texttt{Bridge} reactions. If \texttt{Model} selects a \texttt{Bridge} reaction to execute corresponding to a given \texttt{Polymer}, that \texttt{Polymer} must then choose which \texttt{Polymerase} to move.  To move a \texttt{Polymerase}, the \texttt{Polymer} again takes the scaled movement rates (propensities) and selects a \texttt{Polymerase} randomly, weighted by its propensity. The \texttt{Polymer} attempts to move that \texttt{Polymerase} one position forward. 

\texttt{Polymer} checks for several interactions before \texttt{Polymerase} movement is finalized. First, it checks for a collision with a downstream \texttt{Polymerase} objects by comparing coordinates of the newly-moved \texttt{Polymerase} an upstream \texttt{Polymerase}. If the coordinates overlap, a collision has a occurred. The newly-moved \texttt{Polymerase} moves back to its original starting position, and the current SSA iteration ends. If no collision occurs, \texttt{Polymer} compares the coordinates of the newly-moved \texttt{Polymerase} with that of any upstream \texttt{Terminator}, \texttt{Promoter}, or \texttt{Mask} objects. If the \texttt{Polymerase} overlaps with \texttt{Terminator}, the \texttt{Polymer} verifies that the \texttt{Polymerase} interacts with the \texttt{Terminator}. To simulate readthrough, the \texttt{Polymer} randomly generates a number between 0 and 1. If this value is larger than the the readthrough probability of the \texttt{Terminator}, the \texttt{Polymer} finalizes the movement of \texttt{Polymerase} sets a readthrough flag. This readthrough flag stops \texttt{Polymer} from repeatedly verifying \texttt{Terminator}-\texttt{Polymerase} interactions as the \texttt{Polymerase} moves over the terminator during future SSA iterations. Once \texttt{Polymerase} clears the \texttt{Terminator} completely, \texttt{Polymer} resets the readthrough flag. If no readthrough occurs, \texttt{Polymer} terminates \texttt{Polymerase} by removing the object from the vector of \texttt{Polymerase} objects. The \texttt{Polymer} recomputes its propensity and fires a termination signal to other components of the simulation. This termination signal varies for \texttt{Genome} and \texttt{Transcript} described below.

If \texttt{Polymerase} overlaps with a \texttt{Promoter}, \texttt{Polymer} marks that \texttt{Promoter} as covered and inaccessible. Once \texttt{Polymerase} clears the \texttt{Promoter}, \texttt{Polymer} marks the \texttt{Promoter} as accessible again. \texttt{Polymer} maintains a vector of unbound \texttt{Promoter} objects, and \texttt{SpeciesTracker} maintains a map of which \texttt{Promoter} objects bind to which \texttt{Polymerase} objects. The \texttt{Polymer} reports to \texttt{SpeciesTracker} and \texttt{Model} the number of unbound \texttt{Promoter} objects. If the \texttt{Model} determines that a \texttt{Polymerase} should bind to a \texttt{Promoter}, the \texttt{Polymer} randomly selects the appropriate \texttt{Promoter} to bind, and the polymerase is added to the vector of \texttt{Polymerases} at the \texttt{Promoter} object's position. The newly-bound \texttt{Promoter} is now ready to move on the \texttt{Polymer}.

Lastly, \texttt{Polymerase} objects may interact with \texttt{Mask} objects upon moving. Each \texttt{Polymer} may have a single \texttt{Mask} object. The \texttt{Mask} objects makes portions of the \texttt{Polymer} inaccessible to \texttt{Polymerases}. \texttt{Polymer} treats the \texttt{Mask} as a large \texttt{Polymerase} that may cover the entire \texttt{Polymer}. Upon \texttt{Polymerase} movement, if the \texttt{Polymerase} collides with a \texttt{Mask}, the \texttt{Polymerase} may move back one step, or the \texttt{Mask} may recede. Which of these two interactions occurs depends on the specific \texttt{Mask} and \texttt{Polymerase}, and these interactions differ for \texttt{Genome} and \texttt{Transcript} objects. 

`Genome` and `Transcript`
^^^^^^^^^^^^^^^^^^^^^^^^^

The \texttt{Genome} and \texttt{Transcript} classes are specialized versions of the parent \texttt{Polymer} class (Fig.~\ref{fig:polymer}). A \texttt{Genome} object has a vector member variable that defines a complete transcript template. When a polymerase binds to a promoter, it immediately generates a complete \texttt{Transcript} object based on the transcript template. The newly-generated \texttt{Transcript} object contains genes corresponding to where the polymerase bound and extending to the end of the genome. Upon binding, the polymerase creates a \texttt{Mask} covering the entire \texttt{Transcript}, except for the very 5' end. As the polymerase moves forward from this promoter in the 5'-to-3' direction on the \texttt{Genome} object, the polymerase signals to the \texttt{Transcript} to shift the 5'-end of \texttt{Mask} one base pair towards the 3'-end, thus exposing more of the \texttt{Transcript}. This unmasking process simulates transcript synthesis. Moreover, the length of the transcript, corresponding to the position of the mask, can be determined dynamically as the simulation progresses. The termination position of the transcript does not need to be specified upon promoter binding.

The \texttt{Mask} in a \texttt{Transcript} is inaccessible to ribosomes. Ribosomes, represented by \texttt{Polymerase} objects, collide with the \texttt{Mask} in much the same way that they collide with each other. These collisions simulate ribosomes colliding with an RNA polymerase that is actively synthesizing the transcript on which the ribosome is translating. If a ribosome collides with a mask, the ribosome stalls, just as if it had collided with another ribosome.

A \texttt{Genome} may also define a \texttt{Mask}. This \texttt{Mask} makes portions of the \texttt{Genome} inaccessible to polymerase binding. However, some polymerases are capable of shifting the mask upon colliding with it. This shifting simulates some viral genomes in which a polymerase itself pulls the genome into a cell as it transcribes.

Elements
--------

\texttt{Element} objects are defined as any fixed element along a \texttt{Polymer}. These include promoters, terminators, ribosome binding sites, and stop codons. \texttt{Element} objects may interact with any number of different types of \texttt{Polymerases}. They are also capable of being covered by a \texttt{Polymerase} and thus inaccessible. All \texttt{Element} objects differ from \texttt{Polymerase} and \texttt{Mask} objects in that they have fixed stop and start coordinates. 

Signaling mechanisms
--------------------

A \texttt{Signal} class provides a standardized interface for communication between different objects in Pinetree. For example, when a \texttt{Polymerase} moves it may signal to a transcript \texttt{Mask} that it should also move. When a ribosome reaches a stop codon, it signals to \texttt{SpeciesTracker} that a termination event has occurred and new protein must be added to the species pool. The \texttt{Signal} class follows a ``signals and slots'' model. Some objects carry their own \texttt{Signal} objects. Any function from any object may register with the \texttt{Signal} and occupy one of the \texttt{Signal} slots. These slots represent listeners. An object may then fire a \texttt{Signal} object, transmitting signals to any number of listeners without knowing how such signals will be handled when they reach the listener. This encapsulation allows portions of the simulation to be tested independently from one another. 

Initialization
--------------

A \texttt{Parser} parses the YAML parameter files into Pinetree objects. \texttt{Parser} registers each \texttt{Genome} and \texttt{Reaction} with \texttt{Model} and sometimes \texttt{SpeciesTracker}.

Output
------

One each SSA iteration, \texttt{Model} checks to see if the simulation time has reached the user-specified time step, then it writes information about free species counts from \texttt{SpeciesTracker} to an file.

\begin{figure}[h]
    \centering
    \includegraphics[width=\textwidth]{figures/sim_plots_ms.pdf}
    \caption{Three-gene plasmid simulations demonstrating gene regulation. (A) A single promoter drives transcription of three genes, one of which encodes the RNA polymerase that binds to its own promoter. Differing lengths among genes result in different transcript and protein abundances. No transcripts degrade during the simulation. (B) Two different promoters, responding to different RNA polymerases drive expression of three different genes. The downstream promoter responds to an RNA polymerase encoded upstream (rnapol). This promoter is much stronger than the early promoter, thus producing higher abundances of proteinY transcript and proteinY protein. Both simulations demonstrate non-steady state dynamics.}
    \label{fig:plasmid}
\end{figure}

\begin{figure}[h]
    \centering
    \includegraphics[width=\textwidth]{figures/recoded_plots_ms.pdf}
    \caption{Three-gene plasmid simulations with one gene recoded to use rare codons. (A) Protein abundances of a three-gene model driven by a single promoter, where proteinY has been recoded with rare codons. Translation of proteinY is slower than that of proteinX and rnapol. (B) Ribosome densities on the transcript of proteinY are higher in the recoded gene than in the wildtype. }
    \label{fig:recoded}
\end{figure}

\begin{figure}[h]
    \centering
    \includegraphics[width=\textwidth]{figures/bridge.pdf}
    \caption{Relationship between \texttt{Model}, \texttt{SpeciesTracker}, and individual \texttt{Polymer} objects. Upon execution, the \texttt{Bind} reaction constructs a \texttt{Polymerase} on a \texttt{Polymer} and removes one copy from \texttt{SpeciesTracker}. A \texttt{Bridge} reaction bypasses \texttt{SpeciesTracker} entirely and signals an individual \texttt{Polymer} to move a \texttt{Polymerase}. When a \texttt{Polymerase} reaches a \texttt{Terminator}, the parent \texttt{Polymer} destroys the \texttt{Polymerase} object and increases the copy number count of that \texttt{Polymerase} type in \texttt{SpeciesTracker}. Lastly, \texttt{SpeciesReaction} objects interact only with pooled species in \texttt{SpeciesTracker}.}
    \label{fig:bridge}
\end{figure}

\begin{figure}[h]
    \centering
    \includegraphics[width=0.6\textwidth]{figures/polymer.pdf}
    \caption{Single molecule tracking in Pinetree. When a \texttt{Polymerase} binds to a \texttt{Promoter}, it immediately generates a \texttt{Transcript} object with a \texttt{Mask}. As the \texttt{Polymerase} moves, the \texttt{Mask} object retracts, exposing \texttt{Promoters} on the \texttt{Transcripts}. Dashed lines represent signals between \texttt{Polymerase} and \texttt{Mask} objects. }
    \label{fig:polymer}
\end{figure}