## TCP/IP COMMUNICATION

Table 1 lists the commands expected by the Engine in the TCP/IP communication protocol. There are two types of connections: synchronous (blocked) and asynchronous (non-blocked [NB]) connections. In synchronous connections, the front-end needs to wait until command completion to receive the acknowledge response from the Engine. Asynchronous connections need to be established when the front-end cannot freeze the execution to wait for the acknowledge response, as in virtual reality scenario frontends. This situation happens when a time-consuming command needs to be executed, like "PREPROC" or "TRAIN". In this case, the frontend regularly queries for command termination until the expected response is obtained. To appropriately handle asynchronous connections, a multithread approach with at least two threads is adamant: one, the main thread, that executes all the raw real time processing; and the other, the response thread, which responds to queries of various types of information related to the main thread processing. To make these two threads interoperate, the notion of sessions is introduced. A session is an independent location in the memory of the computer running the engine, capable of storing all the information needed to be sent back to the frontend, such as neurofeedback information and motion corrected volume parameters. The two threads can respond properly to the frontend by a shared access to the same session. The object that handles this TCP/IP communication is the friendEngine object, defined in engine.h file (FES/FES/FRIEND\_Engine/engine.h).


> Table 1. List of commands expected by the Engine in the TCP/IP communication protocol.

| Command | Action |
| --- | --- |
| PREPROC / NBPREPROC | Performs the initial preprocessing steps of FRIEND |
| PIPELINE / NBPIPELINE / FEEDBACK / NBFEEDBACK | Performs the processing of each volume in a run as soon as it becomes available. The difference between the PIPELINE and FEEDBACK commands is that FEEDBACK automatically calculates the feedback values and stores it in a session workspace. |
| GLM / NBGLM | Performs the General Linear Model  (fsl\_glm) calculation of the FSL toolbox calculation |
| FEATURESELECTION / NBFEATURESELECTION | Performs feature selection (see Sato et al., 2013) |
| TRAIN / NBTRAIN | Calls the train function plug-in. |
| TEST | Calls the plug-in feedback calculation function. |
| PLUGIN | Defines the library and the associated plug-in functions to be used in further calls. |
| NEWSESSION | Creates a new session (workspace) in the Engine memory. |
| SESSION/GRAPHPARS | Queries for the movement parameters of a given volume. |
| SESSION/TEST | Returns the feedback information of a volume. |
| SESSION/PREPROC, SESSION/FEEDBACK | Queries if a command, e.g. PREPROC or FEEDBACK has ended. |
| READCONFIG | Sends an entire configuration file associated with the front-end neurofeedback study to set the parameters of the experiment. |