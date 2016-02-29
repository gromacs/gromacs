#include <string>
    
//! Enum to describe the phase corresponding to a property
enum ePhase { epGAS, epLIQUID, epSOLID, epPLASMA, epNR };

/*! \brief
 * Yield string corresponding to phase
 * 
 * \param[in] ep The phase enum
 * \return string corresponding to phase
 */
std::string phase2string(ePhase ep);

/*! \brief
 * Yield phase corresponding to string
 * 
 * \param[in] phase String corresponding to phase
 * \return The phase enum
 */
ePhase string2phase(std::string phase);

