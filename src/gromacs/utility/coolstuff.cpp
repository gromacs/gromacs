/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief Functionality for printing cool strings
 *
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/coolstuff.h"

#include "config.h"

#include <cstdlib>
#include <ctime>

#include <random>
#include <string>

/* This file is completely threadsafe - keep it that way! */

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

//! Whether printing of cool quotes is enabled
bool beCool()
{
    /* Yes, it is bad to check the environment variable every call,
     * but we dont call this routine often, and it avoids using
     * a mutex for locking the variable...
     */
#if GMX_COOL_QUOTES
    return (getenv("GMX_NO_QUOTES") == nullptr);
#else
    /*be uncool*/
    return false;
#endif
}

//! Return a valid random index into \c arrayRef
template<typename T>
const T& getPseudoRandomElement(gmx::ArrayRef<const T> arrayRef)
{
    std::mt19937_64                       generator(std::time(nullptr));
    std::uniform_int_distribution<size_t> distribution(0, arrayRef.size() - 1);
    return arrayRef[distribution(generator)];
}

} // namespace

std::string bromacs()
{
    const char* const bromacsArray[] = {
        "Good gRace! Old Maple Actually Chews Slate",
        "GRoups of Organic Molecules in ACtion for Science",
        "GRowing Old MAkes el Chrono Sweat",
        "Gyas ROwers Mature At Cryogenic Speed",
        "Great Red Owns Many ACres of Sand ",
        "GROningen MAchine for Chemical Simulation",
        "GROup of MAchos and Cynical Suckers",
        "GROtesk MACabre and Sinister",
        "GROwing Monsters And Cloning Shrimps",
        "Great Red Oystrich Makes All Chemists Sane",
        "Good ROcking Metal Altar for Chronical Sinners",
        "Gnomes, ROck Monsters And Chili Sauce",
        "S  C  A  M  O  R  G",
        "Giant Rising Ordinary Mutants for A Clerical Setup",
        "Gromacs Runs On Most of All Computer Systems",
        "Grunge ROck MAChoS",
        "Green Red Orange Magenta Azure Cyan Skyblue",
        "GROningen Mixture of Alchemy and Childrens' Stories",
        "Guyana Rwanda Oman Macau Angola Cameroon Senegal",
        "God Rules Over Mankind, Animals, Cosmos and Such",
        "Georgetown Riga Oslo Madrid Amsterdam Chisinau Stockholm",
        "Gallium Rubidium Oxygen Manganese Argon Carbon Silicon",
        "Glycine aRginine prOline Methionine Alanine Cystine Serine",
        "Gravel Rubs Often Many Awfully Cauterized Sores",
        "Getting the Right Output Means no Artefacts in Calculating Stuff",
        "Gromacs Runs One Microsecond At Cannonball Speeds",
    };

    if (beCool())
    {
        return getPseudoRandomElement<const char*>(bromacsArray);
    }
    else
    {
        return "GROMACS";
    }
}

std::string getCoolQuote()
{
    struct Quote
    {
        const char* text;
        const char* author;
    };

    const Quote quoteArray[] = {
        { "If You Want Something Done You Have to Do It Yourself", "Highlander II" },
        { "I Live the Life They Wish They Did", "Tricky" },
        { "Jesus Built My Hotrod", "Ministry" },
        { "Nurture Another Mind, Before Yours Expires", "Arrested Development" },
        { "Hmm, It *Does* Go Well With the Chicken", "Beastie Boys" },
        { "We Can Dance Like Iggy Pop", "Red Hot Chili Peppers" },
        { "It's So Lonely When You Don't Even Know Yourself", "Red Hot Chili Peppers" },
        { "Do the Dog On the Ground", "Red Hot Chili Peppers" },
        { "Don't Push Me, Cause I'm Close to the Edge", "Tricky" },
        { "Don't Push Me, Cause I'm Close to the Edge", "Grandmaster Flash" },
        { "Bum Stikkie Di Bum Stikkie Di Bum Stikkie Di Bum", "R. Slijngaard" },
        { "She's Not Bad, She's Just Genetically Mean", "Captain Beefheart" },
        { "Being Great is Not So Good", "Red Hot Chili Peppers" },
        { "If Life Seems Jolly Rotten, There's Something You've Forgotten !", "Monty Python" },
        { "Your Proposal is Accepted", "Men In Black" },
        { "Don't Grumble, Give a Whistle !", "Monty Python" },
        { "Stop Drinking My Beer !", "The Amps" },
        { "I Calculate My Birthright", "P.J. Harvey" },
        { "You Should Sleep Late Man, It's Much Easier On Your Constitution", "Beastie Boys" },
        { "You're Insignificant", "Tricky" },
        { "Check Your Output", "P. Ahlstrom" },
        { "What Kind Of Guru are You, Anyway ?", "F. Zappa" },
        { "I Had So Many Problem, and Then I Got Me a Walkman", "F. Black" },
        { "I Caught It In the Face", "P.J. Harvey" },
        { "If You Touch Me, You'll Get Shocked", "Beastie Boys" },
        { "This Puke Stinks Like Beer", "LIVE" },
        { "Art For Arts Sake, Money For Gods Sake", "10 CC" },
        { "A Man Needs a Maid", "N. Young" },
        { "No One Could Foresee the End That Came So Fast", "Slayer" },
        { "Stay Cool, This is a Robbery", "Pulp Fiction" },
        { "With a Little Penknife", "Nick Cave" },
        { "In a Deep Deep Well", "Nick Cave" },
        { "I'm Only Faking When I Get It Right", "Soundgarden" },
        { "Sisters Have Always Fascinated Me", "Speech" },
        { "There's No Room For the Weak", "Joy Division" },
        { "All Work and No Play Makes Jack a Dull Boy", "The Shining" },
        { "They Were So Quiet About It", "Pixies" },
        { "Never Get a Chance to Kick Ass", "The Amps" },
        { "Yeah, a Wuzz, Or a Jerk", "F. Black" },
        { "It's Time to Move On", "F. Black" },
        { "It'll Cure Your Asthma Too !", "F. Zappa" },
        { "Out Of Register Space (Ugh)", "Vi" },
        { "May the Force Be With You", "Star Wars" },
        { "You Try to Run the Universe", "Tricky" },
        { "This May Come As a Shock", "F. Black" },
        { "I Wonder, Should I Get Up...", "J. Lennon" },
        { "I Am Testing Your Grey Matter", "Red Hot Chili Peppers" },
        { "Insane In Tha Membrane", "Cypress Hill" },
        { "You Could Make More Money As a Butcher", "F. Zappa" },
        { "I'll Master Your Language, and In the Meantime I'll Create My Own", "Tricky" },
        { "The Stingrays Must Be Fat This Year", "Red Hot Chili Peppers" },
        { "I'm a Wishbone and I'm Breaking", "Pixies" },
        { "You Think That You're Bigger When You Fill the Void", "Urban Dance Squad" },
        { "And It Goes a Little Something Like This", "Tag Team" },
        { "Kissing You is Like Kissing Gravel", "Throwing Muses" },
        { "You Look Better Upside Down", "Throwing Muses" },
        { "Lunatics On Pogo Sticks", "Red Hot Chili Peppers" },
        { "I Could Take You Home and Abuse You", "Magnapop" },
        { "Move Over Hogey Bear", "Urban Dance Squad" },
        { "You Leave Me Dry", "P.J. Harvey" },
        { "Would You Like to Be the Monster Tonight ?", "Captain Beefheart" },
        { "Meet Me At the Coffee Shop", "Red Hot Chili Peppers" },
        { "She Says She Can't Go Home Without a Chaperone", "E. Costello" },
        { "Keep Your Shoes and Socks On, People", "F. Zappa" },
        { "What If None Of Your Dreams Come True ?", "E. Costello" },
        { "Give a Man a Fish", "Arrested Development" },
        { "The Wheels On the Bus Go Round and Round", "J. Richman" },
        { "I Want to Know Right Now", "Meatloaf" },
        { "What's Your Definition Of Dirty ?", "G. Michael" },
        { "Here's the Way It Might End", "G. Michael" },
        { "Breaking the Law, Breaking the Law", "Judas Priest" },
        { "Just Because the Sun Wants a Place In the Sky", "F. Zappa" },
        { "Baseball Heroes Only", "P.J. Harvey" },
        { "One Cross Each", "Monty Python" },
        { "I Snipe Like Wesley", "Urban Dance Squad" },
        { "Hold On Like Cliffhanger", "Urban Dance Squad" },
        { "It Costs Too Much If It Costs a Lot", "Magnapop" },
        { "Every Sperm is Sacred", "Monty Python" },
        { "Everybody Lie Down On the Floor and Keep Calm", "KLF" },
        { "Love is Like Moby Dick, Get Chewed and Get Spat Out", "Urban Dance Squad" },
        { "Don't Follow Me Home", "Throwing Muses" },
        { "All Beauty Must Die", "Nick Cave" },
        { "I Don't Want to Calm Down", "Throwing Muses" },
        { "We're Gonna Hit You Harder", "Scoter" },
        { "Shake Barrels Of Whisky Down My Throat", "Throwing Muses" },
        { "It's Because Of the Metric System", "Pulp Fiction" },
        { "I Don't Want to Catch Anybody Not Drinking.", "Monty Python" },
        { "This Doesn't Suck, It's a Black Hole !", "K.A. Feenstra" },
        { "Let Me Do This", "Urban Dance Squad" },
        { "I Can't Shake It", "Dinosaur Jr" },
        { "Once Again Let Me Do This", "Urban Dance Squad" },
        { "Pretend That You're Hercule Poirot", "TeX" },
        { "Exactly", "Pulp Fiction" },
        { "Sort Of", "Urban Dance Squad" },
        { "Proceed, With Fingers Crossed", "TeX" },
        { "The Universe is Somewhere In Here", "J.G.E.M. Fraaije" },
        { "You're About to Hurt Somebody", "Jazzy Jeff" },
        { "I Should Be the Pimp Of the Year", "Urban Dance Squad" },
        { "Jesus Can't Save You, Though It's Nice to Think He Tried", "Black Crowes" },
        { "My Heart is Just a Muscle In a Cavity", "F. Black" },
        { "Step Aside, Butch", "Pulp Fiction" },
        { "The World is a Friendly Place", "Magnapop" },
        { "Sometimes Life is Obscene", "Black Crowes" },
        { "Take Your Medications and Preparations and Ram It Up Your Snout", "F. Zappa" },
        { "Count the Bubbles In Your Hair", "The Breeders" },
        { "You Own the Sun", "Throwing Muses" },
        { "I Need a Little Poison", "Throwing Muses" },
        { "Ease Myself Into the Body Bag", "P.J. Harvey" },
        { "Correctomundo", "Pulp Fiction" },
        { "I Don't Like Dirt", "The Breeders" },
        { "Bring Out the Gimp", "Pulp Fiction" },
        { "You Could Be a Shadow", "The Breeders" },
        { "If You're So Special Why aren't You Dead ?", "The Breeders" },
        { "The Path Of the Righteous Man is Beset On All Sides With the Iniquities Of the Selfish "
          "and the Tyranny Of Evil Men.",
          "Pulp Fiction" },
        { "Blessed is He Who In the Name Of Charity and Good Will Shepherds the Weak Through the "
          "Valley Of Darkness, For He is Truly His Brother's Keeper and the Finder Of Lost "
          "Children.",
          "Pulp Fiction" },
        { "And I Will Strike Down Upon Thee With Great Vengeance and With Furious Anger Those Who "
          "Attempt to Poison and Destroy My Brothers.",
          "Pulp Fiction" },
        { "And You Will Know That My Name is the Lord When I Lay My Vengeance Upon Thee.",
          "Pulp Fiction" },
        { "Step On the Brakes", "2 Unlimited" },
        { "You Don't Wanna Know", "Pulp Fiction" },
        { "You Dirty Switch, You're On Again", "The Breeders" },
        { "She's a Good Sheila Bruce", "Monty Python" },
        { "I'm Gonna Get Medieval On Your Ass", "Pulp Fiction" },
        { "Three Little Fonzies", "Pulp Fiction" },
        { "It's Not Your Fault", "Pulp Fiction" },
        { "You Will Be Surprised At What Resides In Your Inside", "Arrested Development" },
        { "The Carpenter Goes Bang Bang", "The Breeders" },
        { "Everybody Wants to Be Naked and Famous", "Tricky" },
        { "Royale With Cheese", "Pulp Fiction" },
        { "Shit Happens", "Pulp Fiction" },
        { "You Fill Your Space So Sweet", "F. Apple" },
        { "Push It Real Good", "Salt 'n' Pepa" },
        { "Check Your Input", "D. Van Der Spoel" },
        { "Catholic School Girls Rule", "Red Hot Chili Peppers" },
        { "It Was My Pleasure", "Pulp Fiction" },
        { "We Don't Bother Anyone", "LIVE" },
        { "I Wrapped a Newspaper Round My Head", "F. Zappa" },
        { "Kick the Dog and You Will Die", "Magnapop" },
        { "We All Get the Flu, We All Get Aids", "LIVE" },
        { "One Ripple At a Time", "Bianca's Smut Shack" },
        { "We Have No Money", "E. Clementi" },
        { "Carry Me Away", "Motors" },
        { "I Solve Problems", "Pulp Fiction" },
        { "A Protein is a Set Of Coordinates", "A.P. Heiner" },
        { "It Doesn't Have to Be Tip Top", "Pulp Fiction" },
        { "Everybody's Good Enough For Some Change", "LIVE" },
        { "It's Against the Rules", "Pulp Fiction" },
        { "I'm An Oakman", "Pulp Fiction" },
        { "I Ripped the Cord Right Out Of the Phone", "Capt. Beefheart" },
        { "I Smell Smoke From a Gun Named Extinction", "Pixies" },
        { "With a Lead Filled Snowshoe", "F. Zappa" },
        { "Right Between the Eyes", "F. Zappa" },
        { "BioBeat is Not Available In Regular Shops", "P.J. Meulenhoff" },
        { "Rub It Right Across Your Eyes", "F. Zappa" },
        { "Shake Yourself", "YES" },
        { "I Am a Wonderful Thing", "Kid Creole" },
        { "Way to Go Dude", "Beavis and Butthead" },
        { "The Microsecond is Within Reach", "P.J. Van Maaren" },
        { "Microsecond Here I Come", "P.J. Van Maaren" },
        { "Confirmed", "Star Trek" },
        { "If You Don't Like Cool Quotes Check Your GMXRC File", "Your Sysadmin" },
        { "When It Starts to Start It'll Never Stop", "Magnapop" },
        { "I'm a Jerk", "F. Black" },
        { "It Wouldn't Hurt to Wipe Once In a While", "Beavis and Butthead" },
        { "Welcome to the Power Age", "2 Unlimited" },
        { "If You See Me Getting High, Knock Me Down", "Red Hot Chili Peppers" },
        { "The Poodle Bites", "F. Zappa" },
        { "The Poodle Chews It", "F. Zappa" },
        { "I Got a Forty Dollar Bill", "F. Zappa" },
        { "We Look Pretty Sharp In These Clothes", "F. Zappa" },
        { "You Got to Relate to It", "A.E. Torda" },
        { "That Was Pretty Cool", "Beavis" },
        { "That Was Really Cool", "Butthead" },
        { "Hang On to Your Ego", "F. Black" },
        { "Pump Up the Volume Along With the Tempo", "Jazzy Jeff" },
        { "Ramones For Ever", "P.J. Van Maaren" },
        { "Have a Nice Day", "R. McDonald" },
        { "Whatever Happened to Pong ?", "F. Black" },
        { "Make the Floor Burn", "2 Unlimited" },
        { "That Was Cool", "Beavis and Butthead" },
        { "These Gromacs Guys Really Rock", "P.J. Meulenhoff" },
        { "You Hear Footsteps Coming From Behind", "Colossal Cave" },
        { "It is Lunchtime", "A.R. Van Buuren" },
        { "You Crashed Into the Swamps", "Silicon Graphics" },
        { "I Am a Poor Lonesome Cowboy", "Lucky Luke" },
        { "Clickety Clickety Click", "System Manager From Hell" },
        { "Been There, Done It", "Beavis and Butthead" },
        { "Load Up Your Rubber Bullets", "10 CC" },
        { "How Do You Like Your Vacation So Far ?", "Speed 2 - Cruise Control" },
        { "It's So Fast It's Slow", "F. Black" },
        { "Ich Bin Ein Berliner", "J.F. Kennedy" },
        { "Take Dehydrated Water On Your Desert Trips", "Space Quest III" },
        { "Your Country Needs YOU", "U.S. Army" },
        { "Don't Eat That Yellow Snow", "F. Zappa" },
        { "I Do It All the Time", "Magnapop" },
        { "Just Give Me a Blip", "F. Black" },
        { "Garbage Collecting...", "GNU Emacs" },
        { "Cut It Deep and Cut It Wide", "The Walkabouts" },
        { "Beat On the Brat With a Baseball Bat", "The Ramones" },
        { "My Head Goes Pop Pop Pop Pop Pop", "F. Black" },
        { "Hangout In the Suburbs If You've Got the Guts", "Urban Dance Squad" },
        { "I Have a Bad Case Of Purple Diarrhea", "Urban Dance Squad" },
        { "It's Bicycle Repair Man !", "Monty Python" },
        { "I've Got Two Turntables and a Microphone", "B. Hansen" },
        { "I Am the Psychotherapist. Please, Describe Your Problems.", "GNU Emacs" },
        { "Watch Out Where the Huskies Go", "F. Zappa" },
        { "I Was Born to Have Adventure", "F. Zappa" },
        { "Is That a Real Poncho ?", "F. Zappa" },
        { "They're Red Hot", "Red Hot Chili Peppers" },
        { "Your Bones Got a Little Machine", "Pixies" },
        { "Oh My God ! It's the Funky Shit", "Beastie Boys" },
        { "Throwing the Baby Away With the SPC", "S. Hayward" },
        { "Engage", "J.L. Picard" },
        { "Everybody is Smashing Things Down", "Offspring" },
        { "Hey Man You Know, I'm Really OK", "Offspring" },
        { "I'm Not Gonna Die Here !", "Sphere" },
        { "I'd Like Monday Mornings Better If They Started Later", "Garfield" },
        { "Here's Another Useful Quote", "S. Boot" },
        { "Wild Pointers Couldn't Drag Me Away", "K.A. Feenstra" },
        { "Let's Go Hang Out In a Mall", "LIVE" },
        { "These are Ideas, They are Not Lies", "Magnapop" },
        { "Bad As This Shit Is, This Shit Ain't As Bad As You Think It Is.", "Jackie Brown" },
        { "My Ass May Be Dumb, But I Ain't No Dumbass.", "Jackie Brown" },
        { "Jesus Not Only Saves, He Also Frequently Makes Backups.", "Myron Bradshaw" },
        { "Player Sleeps With the Fishes", "Ein Bekanntes Spiel Von ID Software" },
        { "Bailed Out Of Edge Synchronization After 10,000 Iterations", "X/Motif" },
        { "God is a DJ", "Faithless" },
        { "Encountered Subspace Anomaly", "Star Trek" },
        { "If I Were You I Would Give Me a Break", "F. Black" },
        { "She Needs Cash to Buy Aspirine For Her Pain", "LIVE" },
        { "Got Coffee, Got Donuts, Got Wasted", "F. Black" },
        { "Boom Boom Boom Boom, I Want You in My Room", "Venga Boys" },
        { "Right Now My Job is Eating These Doughnuts", "Bodycount" },
        { "Wait a Minute, aren't You.... ? (gunshots) Yeah.", "Bodycount" },
        { "If I Wanted You to Understand This, I Would Explain it Better", "J. Cruijff" },
        { "Uh-oh", "Tinky Winky" },
        { "Uh-oh, We're In Trouble", "Shampoo" },
        { "Can't You Make This Thing Go Faster ?", "Black Crowes" },
        { "Get Down In 3D", "George Clinton" },
        { "Uh-oh .... Right Again", "Laurie Anderson" },
        { "(That makes 100 errors; please try again.)", "TeX" },
        { "O My God, They Killed Kenny !", "South Park" },
        { "Drugs are Bad, mmokay", "South Park" },
        { "Let's Unzip And Let's Unfold", "Red Hot Chili Peppers" },
        { "I'd Be Water If I Could", "Red Hot Chili Peppers" },
        { "Space May Be the Final Frontier, But It's Made in a Hollywood Basement",
          "Red Hot Chili Peppers" },
        { "Everything Must Go", "Red Hot Chili Peppers" },
        { "There's Nothing We Can't Fix, 'coz We Can Do It in the Mix", "Indeep" },
        { "It's Coming Right For Us !", "South Park" },
        { "Disturb the Peace of a John Q Citizen", "Urban Dance Squad" },
        { "Wicky-wicky Wa-wild West", "Will Smith" },
        { "This is Tense !", "Star Wars Episode I The Phantom Menace" },
        { "Fly to the Court of England and Unfold", "Macbeth, Act 3, Scene 6, William Shakespeare" },
        { "Why, how now, Claudio ! Whence Comes this Restraint ?",
          "Lucio in Measure for measure, Act 1, Scene 4, William Shakespeare" },
        { "In the End Science Comes Down to Praying", "P. v.d. Berg" },
        { "I'm Looking for a New Simulation", "Stone Temple Pilots" },
        { "I Quit My Job Blowing Leaves", "Beck" },
        { "Live for Liposuction", "Robbie Williams" },
        { "All You Need is Greed", "Aztec Camera" },
        { "You Can Be Too Early, You Can Be Too Late and You Can Be On Time", "J. Cruijff" },
        { "RTFM", "B. Hess" },
        { "Why Do *You* Use Constraints ?", "H.J.C. Berendsen" },
        { "Why Weren't You at My Funeral ?", "G. Groenhof" },
        { "You Can Always Go On Ricky Lake", "Offspring" },
        { "As Always Your Logic Is Impeccable", "Tuvok" },
        { "set: No match.", "tcsh" },
        { "AH ....Satisfaction", "IRIX imapd" },
        { "I Need Love, Not Games", "Iggy Pop & Kate Pierson" },
        { "It's Not Dark Yet, But It's Getting There", "Bob Dylan" },
        { "I Used To Care, But Things Have Changed", "Bob Dylan" },
        { "Working in the Burger Kings, Spitting on your Onion Rings", "Slim Shady" },
        { "Does All This Money Really Have To Go To Charity ?", "Rick" },
        { "Yeah, uh uh, Neil's Head !", "Neil" },
        { "In the Meantime, Take Care of Yourself aaand Eachother", "J. Springer" },
        { "I Feel a Great Disturbance in the Force", "The Emperor Strikes Back" },
        { "Do You Have a Mind of Your Own ?", "Garbage" },
        { "I'll Match Your DNA", "Red Hot Chili Peppers" },
        { "All I Ever Wanted Was Your Life", "Red Hot Chili Peppers" },
        { "Just a Minute While I Reinvent Myself", "Red Hot Chili Peppers" },
        { "There's Still Time to Change the Road You're On", "Led Zeppelin" },
        { "Baby, It Aint Over Till It's Over", "Lenny Kravitz" },
        { "It Just Tastes Better", "Burger King" },
        { "'Nay. We are but men.' Rock!", "Tenacious D" },
        { "Cowardly refusing to create an empty archive", "GNU tar" },
        { "Shaken, not Stirred", "J. Bond" },
        { "Oh, There Goes Gravity", "Eminem" },
        { "Is This the Right Room for an Argument ?", "Monty Python" },
        { "I was detained, I was restrained", "The Smiths" },
        { "The Candlelight Was Just Right", "Beastie Boys" },
        { "Fresh Air, Green Hair", "Frank Black" },
        { "Rat-tat-tat Ka boom boom", "The Smashing Pumpkins" },
        { "Youth is wasted on the young", "The Smashing Pumpkins" },
        { "Miggida-Miggida-Miggida-Mac", "Kriss Kross" },
        { "Interfacing Space and Beyond...", "P. J. Harvey" },
        { "Everything He Lacks, He Makes Up In Denial", "Offspring" },
        { "A Pretty Village Burning Makes a Pretty Fire", "David Sandstrom" },
        { "They don't have any beavers in India, so they have to simulate them", "The Tubes" },
        { "It's Calling Me to Break my Bonds, Again...", "Van der Graaf" },
        { "I believe in miracles cause I'm one", "The Ramones" },
        { "Gabba Gabba Hey!", "The Ramones" },
        { "Shoot them in the back now", "The Ramones" },
        { "Read me your scripture and I will twist it", "Red Hot Chili Peppers" },
        { "Good Music Saves your Soul", "Lemmy" },
        { "Move about like a Scientist, lay down, get kissed", "Red Hot Chili Peppars" },
        { "California, R.I.P.", "Red Hot Chili Peppars" },
        { "Don't You Wish You Never Met Her, Dirty Blue Gene?", "Captain Beefheart" },
        { "Nobody Never Learnt No-Nothing from No History", "Gogol Bordello" },
        { "I'd be Safe and Warm if I was in L.A.", "The Mamas and the Papas" },
        { "It's Unacceptable That Chocolate Makes You Fat", "MI 3" },
        { "My Brothers are Protons (Protons!), My Sisters are Neurons (Neurons)", "Gogol Bordello" },
        { "Put Me Inside SSC, Let's Test Superstring Theory, Oh Yoi Yoi Accelerate the Protons",
          "Gogol Bordello" },
        { "Do You Have Sex Maniacs or Schizophrenics or Astrophysicists in Your Family?",
          "Gogol Bordello" },
        { "Screw a Lightbulb in your Head", "Gogol Bordello" },
        { "Alas, You're Welcome", "Prof. Dumbledore in Potter Puppet Pals" },
        { "Your Shopping Techniques are Amazing", "Gogol Bordello" },
        { "Your Country Raised You, Your Country Fed You, and Just Like Any Other Country it Will "
          "Break You",
          "Gogol Bordello" },
        { "What They Need's a Damn Good Whacking", "The Beatles" },
        { "They Paint Their Faces So Differently From Ours", "Gogol Bordello" },
        { "The Feeling of Power was Intoxicating, Magic", "Frida Hyvonen" },
        { "I was elected to lead, not to read", "President A. Schwarzenegger" },
        { "I managed to get two hours of work done before work", "E. Lindahl" },
        { "Go back to the rock from under which you came", "Fiona Apple" },
        { "It's just the way this stuff is done", "Built to Spill" },
        { "You Fill Me With Inertia", "The Long Blondes" },
        { "I used to be blond and stupid, but now I dyed it black", "Miss Li" },
        { "Aber wenn der Quarterback kommt, um dir die Brille abzunehmen, sag ihm: Danke, die "
          "bleibt wo sie ist",
          "Wir sind Helden" },
        { "Jede der Scherben spiegelt das Licht", "Wir sind Helden" },
        { "Ohne Arbeit waer das Leben oede", "Wir Sind Helden" },
        { "Act like Prometheus would", "Gogol Bordello" },
        { "Making merry out of nothing, like in refugee camp", "Gogol Bordello" },
        { "History has expired", "PubMed Central" },
        { "There's only music to make new ringtones", "Arctic Monkeys" },
        { "Can someone please tell Icarus that he's not the only one falling from the sky?",
          "Urban Dance Squad" },
        { "Ich war schwanger, mir gings zum kotzen", "Nina Hagen" },
        { "What if you're wrong about the great Ju Ju at the bottom of the sea?",
          "Richard Dawkins" },
        { "Come on boys, Let's push it hard", "P.J. Harvey" },
        { "Look at these, my work-strong arms", "P.J. Harvey" },
        { "Is it the invisible chemistry stuff?", "Frida Hyvonen" },
        { "Nada e organico, e tudo programado", "Pitty" },
        { "Sitting on a rooftop watching molecules collide", "A Camp" },
        { "Though the path of the comet is sure, it's constitution is not", "Peter Hammill" },
        { "Everything's formed from particles", "Van der Graaf Generator" },
        { "The time for theory is over", "J. Hajdu" },
        { "What's the point, yo, what's the spread?", "Red Hot Chili Peppers" },
        { "If There Is No Guitar In The House, You Know It's Owner Can Not Be Trusted",
          "Gogol Bordello" },
        { "Carbohydrates is all they groove", "Frank Zappa" },
        { "Never, I said never, compare with experiment", "Magnus Bergh" },
        { "Suzy is a headbanger, her mother is a geek", "The Ramones" },
        { "Now it's filled with hundreds and hundreds of chemicals", "Midlake" },
        { "If it weren't for bad luck, we'd have no luck at all", "The Unthanks" },
        { "There's no way you can rely on an experiment", "Gerrit Groenhof" },
        { "I like to wait, then I feel like I do something", "Carl Caleman" },
        { "Can I have everything louder than everything else?", "Deep Purple" },
        { "He's using code that only you and I know", "Kate Bush" },
        { "Chemical gases filling lungs of little ones", "Black Eyed Peas" },
        { "I've basically become a vegetarian since the only meat I'm eating is from animals I've "
          "killed myself",
          "Mark Zuckerberg" },
        { "Years of calculations and the stress, My science is waiting, nearly complete",
          "Midlake" },
        { "error: too many template-parameter-lists", "g++" },
        { "Science Won't Change You", "The Talking Heads" },
        { "It Doesn't Seem Right, No Computers in Sight", "Faun Fables" },
        { "Some People Say Not to Worry About the Air", "The Talking Heads" },
        { "It seemed a good idea at first", "Gerrit Groenhof" },
        { "There's no kill like overkill, right?", "Erik Lindahl" },
        { "I removed all the lambda defaults so that users have to think!", "Berk Hess" },
        { "I originally implemented PME to prove that you didn't need it...", "Erik Lindahl" },
        { "Take what you want, but just what you need for survival", "Joe Jackson" },
        { "When the universe has expanded, time will contract", "Franz Ferdinand" },
        { "This really is a pretty scene, could you ask your kid to smile please?", "Joe Jackson" },
        { "England's dancing days are done", "P. J. Harvey" },
        { "The future still looks good, and you've got time to rectify all the things that you "
          "should",
          "G. Harrison" },
        { "If humanity has fled shivering from the starry spaces, it has become minutely at home "
          "in the interstices of the speck that it inhabits for an instant",
          "George H. Mead" },
        { "The scientific method is an integral part of human intelligence, and when it has once "
          "been set at work it can only be dismissed by dismissing the intelligence itself",
          "George H. Mead" },
        { "Der Ball ist rund, das Spiel dauert 90 minuten, alles andere ist Theorie", "Lola rennt" },
        { "Life in the streets is not easy", "Marky Mark" },
        { "How will I know it's working right?", "MGMT" },
        { "There was no preconception on what to do", "Daft Punk" },
        { "It takes money to make money, they say", "Lou Reed" },
        { "The future always gets twisted and turned", "Lisa o Piu" },
        { "Do not go where the path may lead, go instead where there is no path and leave a trail",
          "Ralph Waldo Emerson" },
        { "I went to Venice and looked at the paintings of Canaletto to understand how he "
          "presented perspective, and it turned out it was an exponential law. If I had published "
          "this, maybe there would be a Karplus law in art theory as well as the Karplus equation "
          "in NMR",
          "Martin Karplus, Nobel lecture 2013" },
        { "Theoretical chemistry has of course always been important and useful ... at least to "
          "theoretical chemists",
          "Sven Lidin" },
        { "I do not believe continuum electrostatics", "Arieh Warshel, Nobel lecture 2013" },
        { "During my undergraduate work I concluded that electrostatics is unlikely to be "
          "important [for enzymes]",
          "Arieh Warshel, Nobel lecture 2013" },
        { "Martin [Karplus] had a green laser, Arieh [Warshel] had a red laser, I have a *blue* "
          "laser",
          "Michael Levitt, Nobel lecture 2013" },
        { "There's so many shades of black", "The Raconteurs" },
        { "Let us not get carried away with our ideas and take our models too seriously",
          "Nancy Swanson" },
        { "Unfortunately, \"simulation\" has become increasingly misused to mean nothing more than "
          "\"calculation\"",
          "Bill Jorgensen" },
        { "Physics is a few rules, and with some handwaving you can make up the rest",
          "Michael Levitt" },
        { "It doesn't pay to make predictions", "Crowded House" },
        { "Strength is just an accident arising from the weakness of others", "Joseph Conrad" },
        { "On the East coast, a purple patch, to show where the jolly pioneers of progress drink "
          "the jolly lager-beer",
          "Joseph Conrad" },
        { "Restraint! What possible restraint?", "Joseph Conrad" },
        { "It was something to at least have a choice of nightmares", "Joseph Conrad" },
        { "You fight, work, sweat, nearly kill yourself, sometimes you do kill yourself, trying to "
          "accomplish something - and you can't.",
          "Joseph Conrad" },
        { "And after some more talk we agreed that the wisdom of rats had been grossly overrated, "
          "being in fact no greater than that of men",
          "Joseph Conrad" },
        { "It's an easy game, just don't let the ball past!", "Szilard Pall" },
        { "The soul? There's nothing but chemistry here", "Breaking Bad" },
        { "You got one part of that wrong. This is not meth.", "Breaking Bad" },
        { "It's easy to remember: a half a kT is equal to five fourths of a kJ/mol.",
          "Anders Gabrielsson" },
        { "Ubiquitin's just a rock", "Berk Hess" },
        { "... an excellent man, almost worthy of such a wife ...",
          "Jane Eyre in Jane Eyre by Charlotte Bronte" },
        { "Humbug! Most things free-born will submit to anything for a salary",
          "Mr. Rochester in Jane Eyre by Charlotte Bronte" },
        { "Like other defaulters, I like to lay half the blame on ill-fortune and adverse "
          "circumstances",
          "Mr. Rochester in Jane Eyre by Charlotte Bronte" },
        { "Either you will be dashed to atoms on crag points, or lifted up and borne by some "
          "master-wave into a calmer current",
          "Charlotte Bronte" },
        { "I ought to warn you, I have no faith", "Jane Eyre in Jane Eyre by Charlotte Bronte" },
        { "... yet the [economic] profession continued to churn out purely theoretical results "
          "without even knowing what facts needed to be explained.",
          "Thomas Piketty" },
        { "Scientists think they are born with logic; God forbid they should study this discipline "
          "with a history of more than two and a half millennia.",
          "Roald Hoffmann" },
        { "In the processing of models we must be especially cautious of the human weakness to "
          "think that models can be verified or validated. Especially one's own.",
          "Roald Hoffmann" },
        { "... and that dream of dreams, a computational model that predicts everything "
          "accurately.",
          "Roald Hoffmann" },
        { "You see it through a charmed medium: you can not discern that the gilding is slime and "
          "the silk draperies cobwebs; that the marble is sordid slate, and the polished woods "
          "mere refuse chips and scale bark.",
          "Mr. Rochester in Jane Eyre by Charlotte Bronte" },
        { "I know poetry is not dead, nor genius lost; nor has Mammon gained power over either, to "
          "bind or slay; they will both assert their existence, their presence, their liberty and "
          "strength again one day.",
          "Jane Eyre in Jane Eyre by Charlotte Bronte" },
        { "Parallel programming is not about elegance!", "Bill Gropp" },
        { "In a talk you have a choice: You can make one point or no points.", "Paul Sigler" },
        { "Where all think alike, no one thinks very much.", "Walter Lippmann" },
        { "The scientist is not the person who always gives the right answers, he is the one who "
          "asks the right questions.",
          "Claude Levi-Strauss" },
        { "A curious aspect of the theory of evolution is that everybody thinks he understands it.",
          "Jacques Monod" },
        { "When a distinguished but elderly scientist states that something is possible, he is "
          "almost certainly right. When he states that something is impossible, he is very "
          "probably wrong.",
          "Arthur C. Clarke" },
        { "Energy is a very subtle concept. It is very, very difficult to get right.",
          "Richard Feynman" },
        { "The determined Real Programmer can write FORTRAN programs in any language.", "Ed Post" },
        { "FORTRAN was the language of choice for the same reason that three-legged races are "
          "popular.",
          "Ken Thompson" },
        { "A computer without COBOL and FORTRAN is like a piece of chocolate cake without ketchup "
          "or mustard.",
          "Unix fortune program" },
        { "Consistently separating words by spaces became a general custom about the tenth century "
          "A.D., and lasted until about 1957, when FORTRAN abandoned the practice.",
          "Sun FORTRAN Reference Manual" },
        { "Ludwig Boltzmann, who spent much of his life studying statistical mechanics, died in "
          "1906, by his own hand. Paul Ehrenfest, carrying on the same work, died similarly in "
          "1933. Now it is our turn to study statistical mechanics. Perhaps it will be wise to "
          "approach the subject cautiously.",
          "David Goodstein" },
        { "It all works because Avogadro's number is closer to infinity than to 10.",
          "Ralph Baierlein" },
        { "In this house, we OBEY the laws of thermodynamics!", "Homer Simpson" },
        { "We mathematicians are all a bit crazy.", "Lev Landau" },
        { "There is no such thing as free energy. Anyone who advocates it does not know what he is "
          "talking about.",
          "Alireza Haghighat" },
        { "In science it often happens that scientists say, 'You know that's a really good "
          "argument; my position is mistaken,' and then they would actually change their minds and "
          "you never hear that old view from them again. They really do it. It doesn't happen as "
          "often as it should, because scientists are human and change is sometimes painful. But "
          "it happens every day. I cannot recall the last time something like that happened in "
          "politics or religion.",
          "Carl Sagan" },
        { "There is nothing new to be discovered in physics now. All that remains is more and more "
          "precise measurement.",
          "Lord Kelvin, 1900" },
        { "I love fools' experiments. I am always making them.", "Charles Darwin" },
        { "If you want to save your child from polio, you can pray or you can inoculate... choose "
          "science.",
          "Carl Sagan" },
        { "Molecular biology is essentially the practice of biochemistry without a license.",
          "Edwin Chargaff" },
        { "If at one time or another I have brushed a few colleagues the wrong way, I must "
          "apologize: I had not realized that they were covered with fur.",
          "Edwin Chargaff" },
        { "It has not escaped our notice that the specific pairing we have postulated immediately "
          "suggests a possible copying mechanism for the genetic material.",
          "Watson & Crick" },
        { "The researcher's art is first of all to find himself a good boss.", "Andre Lwoff" },
        { "What about my nose?",
          "Aneesur Rahman, responding to an Argonne manager arguing the long hair of Charles "
          "Bennett in his group was disreputing the lab; Retold by Michael Klein" },
        { "Science, my lad, is made up of mistakes, but they are mistakes which it is useful to "
          "make, because they lead little by little to the truth.",
          "Jules Verne" },
        { "Don't be afraid of hard work. Nothing worthwhile comes easily. Don't let others "
          "discourage you or tell you that you can't do it. In my day I was told women didn't go "
          "into chemistry. I saw no reason why we couldn't.",
          "Gertrude Elion" },
        { "The Nobel Prize is fine, but the drugs I've developed are rewards in themselves.",
          "Gertrude Elion" },
        { "...sometimes a scream is better than a thesis.", "Ralph Waldo Emerson" },
        { "The great tragedy of science - the slaying of a beautiful hypothesis by an ugly fact.",
          "Thomas Henry Huxley" },
        { "Dr Pauling, how do you have so many good ideas? Well David, I have a lot of ideas and "
          "throw away the bad ones.",
          "Linus Pauling" },
        { "I try to identify myself with the atoms... I ask what I would do If I were a carbon "
          "atom or a sodium atom.",
          "Linus Pauling" },
        { "I admired Bohr very much. We had long talks together, long talks in which Bohr did "
          "practically all the talking.",
          "Paul Dirac" },
        { "Predictions can be very difficult - especially about the future.", "Niels Bohr" },
        { "For those who want some proof that physicists are human, the proof is in the idiocy of "
          "all the different units which they use for measuring energy.",
          "Richard Feynman" },
        { "Dreams seldom materialize on their own.", "Dian Fossey" },
        { "Above all, don't fear difficult moments. The best comes from them.",
          "Rita Levi-Montalcini" },
        { "Our struggle today is not to have a female Einstein get appointed as an assistant "
          "professor. It is for a woman schlemiel to get as quickly promoted as a male schlemiel.",
          "Bella Abzug" },
        { "I never thought of stopping, and I just hated sleeping. I can't imagine having a better "
          "life.",
          "Barbara McClintock" },
        { "The farther the experiment is from theory, the closer it is to the Nobel Prize.",
          "Irene Joliot-Curie" },
        { "I never see what has been done; I only see what remains to be done.", "Marie Curie" },
        { "There is no reason for any individual to have a computer in his home.",
          "Ken Olsen, head of Digital Equipment Corp." },
        { "People disagree with me. I just ignore them.",
          "Linus Torvalds on the use of C++ in the kernel" },
        { "Beware of bugs in the above code; I have only proved it correct, not tried it.",
          "Donald Knuth" },
        { "My greatest contribution to the field of science is that I never entered it.",
          "Colin Powell" },
        { "We are perhaps not far removed from the time when we shall be able to submit the bulk "
          "of chemical phenomena to calculation.",
          "Joseph Gay-Lussac, 1808" },
        { "If mathematical analysis should ever hold a prominent place in chemistry - an "
          "aberration which is happily almost impossible - it would occasion a rapid and "
          "widespread degeneration of that science.",
          "Aguste Comte, 1830" },
        { "Almost without exception, the talented women I have known have believed they had less "
          "ability than they actually had. And almost without exception, the talented men I have "
          "known believed they had more.",
          "Gregory Petsko" },
        { "The first 90% of the code accounts for the first 90% of the development time. The "
          "remaining 10% of the code accounts for the other 90% of the development time.",
          "Tom Cargill" },
        { "The Internet?  We are not interested in it.", "Bill Gates, 1993" },
        { "Perl: The only language that looks the same before and after RSA encryption.",
          "Keith Bostic" },
        { "There are only two things wrong with C++:  The initial concept and the implementation.",
          "Bertrand Meyer" },
        { "XML is not a language in the sense of a programming language any more than sketches on "
          "a napkin are a language.",
          "Charles Simonyi" },
        { "It has been discovered that C++ provides a remarkable facility for concealing the "
          "trivial details of a program - such as where its bugs are.",
          "David Keppel" },
        { "UNIX is basically a simple operating system. It just takes a genius to understand its "
          "simplicity.",
          "Dennis Ritchie" },
        { "There are only two kinds of programming languages: those people always bitch about and "
          "those nobody uses.",
          "Bjarne Stroustrup" },
        { "If Java had true garbage collection, most programs would delete themselves upon "
          "execution.",
          "Robert Sewell" },
        { "Documentation is like sex: When it's good it's great, and when it's bad it's better "
          "than nothing.",
          "Linus Torvalds" },
        { "C has the power of assembly language and the convenience of... assembly language.",
          "Dennis Ritchie" },
        { "The last good thing written in C was Franz Schubert's Symphony Number 9.",
          "Erwin Dieterich" },
        { "User-friendly, adj.: Programmer-hostile.", "New Hacker's Dictionary" },
        { "First off, I'd suggest printing out a copy of the GNU coding standards, and NOT read "
          "it. Burn them, it's a great symbolic gesture.",
          "Linus Torvalds" },
        { "I invented the term 'Object-Oriented', and I can tell you I did not have C++ in mind.",
          "Alan Kay, author of Smalltalk" },
        { "FORTRAN, the infantile disorder, by now nearly 20 years old, is hopelessly inadequate "
          "for whatever computer application you have in mind today: it is now too clumsy, too "
          "risky, and too expensive to use.",
          "Edsger Dijkstra, 1970" },
        { "Do you know what cations don't like? Dog-ions. Do you know what they like? Pie.",
          "Tom Cheatham" },
        { "The most exciting phrase to hear in science, the one that heralds new discoveries, is "
          "not \"Eureka\" but \"That's funny...\".",
          "Isaac Asimov" },
        { "Those people who think they know everything are a great annoyance to those of us who "
          "do.",
          "Isaac Asimov" },
        { "No great discovery was ever made without a bold guess.", "Marie Curie" },
        { "Chance favors the prepared mind.", "Louis Pasteur" },
        { "I love deadlines. I like the whooshing sound they make as they fly by.",
          "Douglas Adams" },
        { "Good judgement is the result of experience; experience is the result of bad judgement.",
          "Mark Twain" },
        { "No matter how important you are, you are not as important as lunch.", "Randy Pausch" },
        { "There is just one thing I can promise you about the outer-space program: your tax "
          "dollar will go farther.",
          "Wernher von Braun" },
        { "Harvard makes mistakes too, you know. Kissinger taught there.", "Woody Allen" },
        { "Nothing in biology makes sense except in the light of evolution.",
          "Theodosius Dobzhansky" },
        { "I have a hunch that the unknown sequences of DNA will decode into copyright notices and "
          "patent protections.",
          "Donald Knuth" },
        { "It always takes longer than you think even when you take Hofstadter's Law into account.",
          "Hofstadter's Law" },
        { "A ship in port is safe, but that is not what ships are for. Sail out to sea and do new "
          "things.",
          "Grace Hopper, developer of COBOL" },
        { "I was told I'd never make it to VP rank because I was too outspoken. Maybe so, but I "
          "think men will always find an excuse for keeping women in their 'place.' So, let's make "
          "that place the executive suite and start more of our own companies.",
          "Jean Bartik, ENIAC developer" },
        { "If it's a good idea, go ahead and do it. It's much easier to apologize than it is to "
          "get permission.",
          "Grace Hopper, developer of COBOL" },
        { "This isn't right. This isn't even wrong.", "Wolfgang Pauli" },
        { "Louis Pasteur's theory of germs is ridiculous fiction.",
          "Pierre Pachet, Professor of Physiology at Toulouse, 1872" },
        { "Research ! A mere excuse for idleness; it has never achieved, and will never achieve "
          "any results of the slightest value.",
          "Benjamin Jowett, British theologian, 1817-93" },
        { "Problems worthy of attack prove their worth by hitting back.", "Piet Hein" },
        { "You should never bet against anything in science at odds of more than about 10^12 to 1.",
          "Ernest Rutherford" },
        { "X-rays will prove to be a hoax.", "Lord Kelvin, while president of the Royal Society" },
        { "If you're doing I/O, you're doing it wrong!", "Cannada \"Drew\" Lewis" },
        { "The easiest way to scale well is to have bad single-core performance", "Blind Freddie" },
        { "Heard a talk introducing a new language called Swift, from a guy named Wozniak, and it "
          "had nothing to do with Apple!",
          "Adam Cadien" },
        { "When doing HPC, don't communica", "Jim Demmel" },
        { "Today we're not going to optimize our CUDA code, cause that's just a rabbit hole of "
          "misery!",
          "Tim Warburton" },
        { "Big Data is like teenage sex: everyone talks about it, nobody really knows how to do "
          "it, everyone thinks everyone else is doing it, so everyone claims they are doing it...",
          "Dan Ariely" },
        { "It seems likely that significant software contributions to existing scientific software "
          "projects are not likely to be rewarded through the traditional reputation economy of "
          "science. Together these factors provide a reason to expect the over-production of "
          "independent scientific software packages, and the underproduction of collaborative "
          "projects in which later academics build on the work of earlier ones.",
          "Howison & Herbsleb" },
        { "On average, it takes twenty years for the world's largest super computer to shrink down "
          "to the size of your laptop.",
          "Pete Beckman" },
        { "When using an abacus, a human can achieve about 0.1 flops/watt. Super-computers achieve "
          "about 2 gigaflops/watt.",
          "John Linford" },
        { "Try to calculate the numbers that have been", "The Smoke Fairies" },
        { "Please implement proper hep writing", "GROMACS" },
        { "The three principal virtues of a programmer are Laziness, Impatience, and Hubris",
          "Larry Wall" },
        { "You're like them scientists on TV explaining black holes. More you talk, less I get",
          "Jess Walter" },
        { "Wedged as we are between two eternities of idleness, there is no excuse for being idle "
          "now",
          "Anthony Burgess" },
        { "Even the *healthy* people move in clouds of cigarette smoke, women straining polyester, "
          "men in raggedly cutoffs slathering mayonnaise on foot-long hot dogs. It's as if the "
          "hotel were hosting a conference on adult onset diabetes",
          "Jess Walter" },
        { "In practice, throwing traditional norms and values overboard results not in perfect "
          "freedom and relationships based on reason, but in chaos and fear",
          "Paul Verhaeghe" },
        { "When I asked a younger colleague at the university how he had been able to change his "
          "research field several times within a decade or so, he answered: \"It's just a question "
          "of new software\"",
          "Paul Verhaeghe" },
        { "Never mind, death professor, your structure's fine", "TV on the Radio" },
        { "Come and play on the hospital roof, I got something that's yours", "Sherlock" },
        { "Njuta men inte frossa, springa men inte fly", "Paganus" },
        { "Misslycka kan man med all kod", "Mats Nylen" },
        { "Two guys can move very fast when they're motivated enough and unemployed",
          "Eric Betzig" },
        { "A protein is a chain of letters.", "Julie Bernauer" },
        { "The best way to obtain plausible negative examples is to run a docking program with a "
          "biophysics-based function.",
          "Julie Bernauer" },
        { "I think everybody should like everybody.", "Andy Warhol" },
        { "But I always say, one's company, two's a crowd, and three's a party.", "Andy Warhol" },
        { "We'll celebrate a woman for anything, as long as it's not her talent.",
          "Colleen McCullough" },
        { "I believe the big bang of self-driving cars is about to come.",
          "Jen-Hsun Huang, CEO NVIDIA" },
        { "This is where we have been working hard to push down performance.",
          "Szilard Pall, GTC 2015 talk" },
        { "Some of these pro-drug messages come from popular culture.", "John Walters" },
        { "Don't pay any attention to what they write about you. Just measure it in inches.",
          "Andy Warhol" },
        { "Art is what you can get away with.", "Andy Warhol" },
        { "I spent a lot of money on booze, birds and fast cars. The rest I just squandered.",
          "George Best" },
        { "The only greatness for man is immortality.", "James Dean" },
        { "Do not quench your inspiration and your imagination; do not become the slave of your "
          "model.",
          "Vincent Van Gogh" },
        { "You always pass failure on the way to success.", "Mickey Rooney" },
        { "I always seem to get inspiration and renewed vitality by contact with this great novel "
          "land of yours which sticks up out of the Atlantic.",
          "Winston Churchill" },
        { "I am at two with nature.", "Woody Allen" },
        { "I'm no model lady. A model's just an imitation of the real thing.", "Mae West" },
        { "Science is the great antidote to the poison of enthusiasm and superstition.",
          "Adam Smith, Wealth of Nations, 1776" },
        { "Science is a wonderful thing if one does not have to earn one's living at it.",
          "Albert Einstein" },
        { "Science is the record of dead religions.", "Oscar Wilde" },
        { "Physics isn't a religion. If it were, we'd have a much easier time raising money.",
          "Leon Lederman" },
        { "It is now quite lawful for a Catholic woman to avoid pregnancy by a resort to "
          "mathematics, though she is still forbidden to resort to physics and chemistry.",
          "Henry Louis Mencken" },
        { "An expert is a person who has made all the mistakes that can be made in a very narrow "
          "field.",
          "Niels Bohr" },
        { "In my opinion, we don't devote nearly enough scientific research to finding a cure for "
          "jerks.",
          "Bill Watterson" },
        { "Scientists do not join hands every Sunday and sing \"Yes gravity is real! I know "
          "gravity is real! I will have faith! I believe in my heart that what goes up, up, up "
          "must come down, down, down. Amen!\" If they did, we would think they were pretty "
          "insecure about the concept.",
          "Dan Barker" },
        { "Take away paradox from the thinker and you have a professor.", "Soren Kirkegaard" },
        { "Measuring programming progress by lines of code is like measuring aircraft building "
          "progress by weight.",
          "Bill Gates" },
        { "Protons give an atom its identity, electrons its personality.", "Bill Bryson" },
        { "Money won't buy happiness, but it will pay the salaries of a large research staff to "
          "study the problem.",
          "Bill Vaughan" },
        { "Torture numbers, and they'll confess to anything.", "Greg Easterbrook" },
        { "Should we force science down the throats of those that have no taste for it? Is it our "
          "duty to drag them kicking and screaming into the twenty-first century? I am afraid that "
          "it is.",
          "George Porter" },
        { "A computer would deserve to be called intelligent if it could deceive a human into "
          "believing that it was human.",
          "Alan Turing" },
        { "Any one who considers arithmetical methods of producing random digits is, of course, in "
          "a state of sin.",
          "John von Neumann" },
        { "No, no, you're not thinking, you're just being logical.", "Niels Bohr" },
        { "As an adolescent I aspired to lasting fame, I craved factual certainty, and I thirsted "
          "for a meaningful vision of human life -- so I became a scientist. This is like becoming "
          "an archbishop so you can meet girls.",
          "Matt Cartmill" },
        { "Problems worthy / of attack / prove their worth / by hitting back.", "Piet Hein" },
        { "Naive you are if you believe life favours those who aren't naive.", "Piet Hein" },
        { "Never measure the height of a mountain until you have reached the top. Then you will "
          "see how low it was.",
          "Dag Hammarskjold" },
        { "Praise those of your critics for whom nothing is up to standard.", "Dag Hammarskjold" },
        { "Inventions have long since reached their limit, and I see no hope for further "
          "development.",
          "Julius Sextus Frontinus, 1st century A.D." },
        { "Lottery: A tax on people who are bad at math.", "Ambrose Bierce" },
        { "Even if you are on the right track, you will get run over if you just sit there.",
          "Will Rogers" },
        { "Programming today is a race between software engineers striving to build bigger and "
          "better idiot-proof programs, and the universe trying to build bigger and better idiots. "
          "So far, the universe is winning.",
          "Rick Cook" },
        { "There's a limit to how many times you can read how great you are and what an "
          "inspiration you are, but I'm not there yet.",
          "Randy Pausch" },
        { "Throughout my academic career, I'd given some pretty good talks. But being considered "
          "the best speaker in the computer science department is like being known as the tallest "
          "of the Seven Dwarfs.",
          "Randy Pausch" },
        { "If everything seems under control, you're just not going fast enough.",
          "Mario Andretti" },
        { "Sincerity is the key to success. Once you can fake that you've got it made.",
          "Groucho Marx" },
        { "This work contains many things which are new and interesting. Unfortunately, everything "
          "that is new is not interesting, and everything which is interesting, is not new.",
          "Lev Landau" },
        { "Does college pay? They do if you are a good open-field runner.", "Will Rogers" },
        { "Academe, n.: An ancient school where morality and philosophy were taught. Academy, n.: "
          "A modern school where football is taught.",
          "Ambrose Bierce" },
        { "This simulation is not as the former.",
          "Malvolio, Act II, scene V of Shaphespeare's Twelfth Night" },
        { "Here, kitty, kitty...", "Erwin Schroedinger" },
        { "Sir, spare your threats: The bug which you would fright me with I seek.",
          "Hermione, Act III, scene II of Shakespeare's Winter's Tale" },
        { "Erwin with his psi can do / Calculations quite a few. / But one thing has not been seen "
          "/ Just what psi really mean.",
          "Felix Bloch" },
        { "Only entropy comes easy.", "Anton Chekov" },
        { "The loveliest theories are being overthrown by these damned experiments; it is no fun "
          "being a chemist any more.",
          "Justus von Liebig, letter to J.J. Berzelius 1834" },
        { "If all else fails, immortality can always be assured by spectacular error.",
          "John Kenneth Galbraith" },
        { "Always code as if the person who ends up maintaining your code is a violent psychopath "
          "who knows where you live.",
          "Martin Golding" },
        { "If I have not seen as far as others, it is because giants were standing on my "
          "shoulders.",
          "Hal Abelson" },
        { "Weaseling out of things is important to learn. It's what separates us from the "
          "animals... except the weasels.",
          "Homer Simpson" },
        { "In science, truth always wins.", "Max Perutz" },
        { "Creativity in science, as in art, cannot be organized. It arises spontaneously from "
          "individual talent. Well-run laboratories can foster it, but hierarchical organizations, "
          "inflexible bureaucratic rules, and mountains of futile paperwork can kill it.",
          "Max Perutz" },
        { "Every electron is sacred.", "Greg McMullan, on Cryo-EM detectors" },
        { "Science adjusts its views based on what's observed. Faith is the denial of observation "
          "so that belief can be preserved.",
          "Tim Minchin" },
        { "Isn't this enough? Just this world? Just this beautiful, complex wonderfully "
          "unfathomable world? How does it so fail to hold our attention that we have to diminish "
          "it with the invention of cheap, man-made myths and monsters?",
          "Tim Minchin" },
        { "If you open your mind too much, your brains will fall out.", "Tim Minchin" },
        { "\"Everything organic and natural is good\" - ignoring the fact that organic natural "
          "substances include arsenic and poo and crocodiles. And everything chemical is bad, "
          "ignoring the fact that... everything is chemicals.",
          "Tim Minchin" },
        { "A program that has not been tested does not work.", "Bjarne Stroustrup" },
        { "You could give Aristotle a tutorial. And you could thrill him to the core of his being. "
          "Such is the privilege of living after Newton, Darwin, Einstein, Planck, Watson, Crick "
          "and their colleagues.",
          "Richard Dawkins" },
        { "A robot will be truly autonomous when you instruct it to go to work and it decides to "
          "go to the beach instead.",
          "Brad Templeton" },
        { "If you want to destroy my sweater, hold this thread as I walk away.", "Weezer" },
        { "To survive science you have to become science.", "Gerrit Groenhof" },
        { "Contemplating answers that could break my bonds.", "Peter Hammill" },
        { "I always think there is something foreign about jolly phrases at breakfast.",
          "Mr. Carson in Downtown Abbey" },
        { "According to my computations we're overdue for a transformation.", "Jackson Browne" },
        { "Therefore, things must be learned only to be unlearned again or, more likely, to be "
          "corrected.",
          "Richard Feynman" },
        { "You wouldn't walk into a chemistry lab and mix two clear liquids together just because "
          "they look pretty much the same, would you?",
          "Justin Lemkul" },
        { "They don't have half hours in the north", "Carl Caleman" },
        { "Safety lights are for dudes", "Ghostbusters 2016" },
        { "It's 2040 now. Our President is a plant.", "Ghostbusters 2016" },
        { "It's just B I O L O G Y, can't you see?", "Joe Jackson" },
        { "Input, output, electricity", "Joni Mitchell" },
        { "Your daddy ain't your daddy but your daddy don't know", "Dalahan" },
        { "Why is the Earth moving 'round the sun? Floating in the vacuum with no purpose, not a "
          "one",
          "Fleet Foxes" },
        { "Everybody has a plan until they get punched in the mouth", "Mike Tyson" },
        { "Sacrifices must be made",
          "Otto Lilienthal, dying after having crashed with his glider in 1896" },
        { "The secret to getting ahead is getting started", "Mark Twain" },
        { "Water is just water", "Berk Hess" },
        { "GROMACS First : Making MD Great Again", "Vedran Miletic" },
        { "You still have to climb to the shoulders of the giants", "Vedran Miletic" },
        { "The road to openness is paved with git commits", "Vedran Miletic" },
        { "Performance and power are great targets for tuning, but really you want to tune for "
          "money!",
          "Erik Lindahl" },
        { "Here are all the 'gmx' tools... but no gmx writethesis", "Christian Blau" },
        { "The best part of winter in Stockholm is going to Australia", "Mark Abraham" },
        { "If you don't know what you're doing, use a (M)BAR-based method", "Erik Lindahl" },
        { "All models are wrong, but some are useful.", "George Box" },
        { "If your experiment needs a statistician, you need a better experiment.",
          "Ernest Rutherford" },
        { "Facts are stubborn things, but statistics are more pliable.", "Laurence Peter" },
        { "In ancient times they had no statistics so they had to fall back on lies.",
          "Stephen Leacock" },
        { "If at first you don't succeed, try two more times so that your failure is statistically "
          "significant.",
          "Dallas Warren" },
        { "Your theory is crazy, but it's not crazy enough to be true.", "Niels Bohr" },
        { "Science may never come up with a better office communication system than the coffee "
          "break.",
          "Earl Wilson" },
        { "A scientific truth does not triumph by convincing its opponents and making them see the "
          "light, but rather because its opponents eventually die and a new generation grows up "
          "that is familiar with it.",
          "Max Planck" },
        { "Computer science is no more about computers than astronomy is about telescopes",
          "Edsger Dijkstra" },
        { "If we knew what it was we were doing, it would not be called research, would it?",
          "Albert Einstein" },
        { "I have not failed. I've just found 10,000 ways that won't work", "Thomas Alva Edison" },
        { "The public have an insatiable curiosity to know everything, except what is worth "
          "knowing.",
          "Oscar Wilde" },
        { "Philosophy of science is about as useful to scientists as ornithology is to birds.",
          "Richard Feynman" },
        { "I had trouble with physics in college. When I signed up I thought it said psychics.",
          "Greg Tamblyn" },
        { "There's an old saying among scientific guys: “You can't make an omelet without breaking "
          "eggs, ideally by dropping a cement truck on them from a crane.",
          "Dave Barry" },
        { "Occams Razor is the scientific principle that, all things being equal, the simplest "
          "explanation is always the dog ate my homework.",
          "Greg Tamblyn" },
        { "When you get right down to it, almost every explanation Man came up with for anything "
          "until about 1926 was stupid.",
          "Dave Barry" },
        { "We all understand the twinge of discomfort at the thought that we share a common "
          "ancestor with the apes. No one can embarrass you like a relative.",
          "Neal DeGrasse Tyson" },
        { "In physics, you don't have to go around making trouble for yourself. Nature does it for "
          "you.",
          "Frank Wilczek" },
        { "Every revolutionary idea seems to evoke three stages of reaction. They may be summed up "
          "by the phrases: (1) It's completely impossible. (2) It's possible, but not worth doing. "
          "(3) I said it was a good idea all along.",
          "Arthur C. Clarke" },
        { "Computers are like humans - they do everything except think.", "John von Neumann" },
        { "With four parameters I can fit an elephant, and with five I can make him wiggle his "
          "trunk.",
          "John von Neumann" },
        { "Christianity may be OK between consenting adults in private but should not be taught to "
          "young children.",
          "Francis Crick" },
        { "All approaches at a higher level are suspect until confirmed at the molecular level.",
          "Francis Crick" },
        { "We haven't the money, so we've got to think.", "Ernest Rutherford" },
        { "Furious activity is no substitute for understanding.", "H.H. Williams" },
        { "Discovery: A couple of months in the laboratory can frequently save a couple of hours "
          "in the library.",
          "Anonymous" },
        { "Never replicate a successful experiment.", "Fett's law." },
        { "Raw data is like raw sewage, it requires some processing before it can be spread "
          "around. The opposite is true of theories.",
          "Jim Carr" },
        { "A university faculty is 500 egotists with a common parking problem.", "Keith Sullivan" },
        { "Studying expands knowledge. Knowledge is power. Power corrupts. Corruption is a crime. "
          "Crime doesn't pay.",
          "Anonymous" },
        { "A professor is one who talks in someone else's sleep.", "W.H. Auden" },
        { "A tidy laboratory means a lazy chemist.", "J.J. Berzelius" },
        { "Microbiology Lab - Staph Only.", "Anonymous" },
        { "I can't go to a restaurant and order food because I keep looking at the fonts on the "
          "menu. Five minutes later I realize that it's also talking about food.",
          "Donald Knuth" },
        { "Physics is like sex: sure, it may give some practical results, but that's not why we do "
          "it",
          "Richard P. Feynman" },
        { "Statistics: The only science that enables different experts using the same figures to "
          "draw different conclusions.",
          "Evan Esar" },
        { "If I could remember the names of all these particles, I'd be a botanist.",
          "Albert Einstein" },
        { "Science... never solves a problem without creating ten more.", "George Bernard Shaw" },
        { "A mathematician is a blind man in a dark room looking for a black cat which isn't "
          "there.",
          "Charles Darwin" },
        { "Nothing shocks me. I'm a scientist.", "Harrison Ford as Indiana Jones" },
        { "There is an infinite set A that is not too big.", "John von Neumann" },
        { "If it's all right with Dirac, it's all right with me.",
          "Enrico Fermi, on being told that there was experimental evidence He-3 nuclei obey "
          "Fermi-Dirac statistics." },
        { "I cannot think of a single one, not even intelligence.",
          "Enrico Fermi, when asked what characteristics physics Nobel laureates had in common." },
        { "Heavier-than-air flying machines are impossible.",
          "Lord Kelvin, President of Royal Society, 1895." },
        { "All that glitters may not be gold, but at least it contains free electrons.",
          "John Desmond Baernal" },
        { "It is disconcerting to reflect on the number of students we have flunked in chemistry "
          "for not knowing what we later found to be untrue.",
          "Robert L. Weber" },
        { "People are DNA's way of making more DNA.", "Edward O. Wilson" },
        { "The best model of a cat is another cat..., specially the same cat.",
          "Arturo Rosenblueth" },
        { "Computer dating is fine, if you are a computer.", "Rita May Brown" },
        { "The most likely way for the world to be destroyed, most experts agree, is by accident. "
          "That's where we come in; we're computer professionals. We cause accidents.",
          "Nathaniel Borenstein" },
        { "An intellectual is someone who has found something more interesting than sex.",
          "Edgar Wallace" },
        { "Base eight is just like base ten really, if you're missing two fingers.", "Tom Lehrer" },
        { "If 10 years from now, when you are doing something quick and dirty, you suddenly "
          "visualize that I am looking over your shoulders and say to yourself: 'Dijkstra would "
          "not have liked this', well that would be enough immortality for me.",
          "Edsger Dijkstra" },
        { "Memory is like an orgasm. It's a lot better if you don't have to fake it.",
          "Seymour Cray, on virtual memory" },
        { "A computer once beat me at chess, but it was no match for me at kick boxing.",
          "Emo Philips" },
        { "Home computers are being called upon to perform many new functions, including the "
          "consumption of homework formerly eaten by the dog.",
          "Doug Larson" },
        { "Forcefields are like dating; things go fine for a while and then sometimes it goes "
          "really bad.",
          "Alex MacKerell" },
        { "This type of advanced sampling techniques... which are not so advanced, really.",
          "Viveca Lindahl, on AWH, at her thesis defense." },
        { "C++ is tricky. You can do everything. You can even make every mistake.",
          "Nicolai Josuttis, CppCon2017" },
        { "Why would the backup server database get corrupted anyway?",
          "Stefan Fleischmann -- system administrator, physicist, optimist." },
        { "Teaching quantum computing is like teaching computer science at Hogwarts.",
          "Thomas Sterling, ISC2018 keynote" },
        { "It is unfortunate that the authors did not make better use of all the electric power "
          "energy that went into these massive computations.",
          "An anonymous referee" },
        { "Doctor, doctor, it hurts when I hit myself in the head with the hammer! - So don't do "
          "it!",
          "Bjarne Stroustrup at CppCon2015" },
        { "This is extremely unlikely.", "Berk Hess" },
        { "Nothing is more anarchic than power.", "Pier Paolo Pasolini" },
        { "Never attribute to malice that which can be adequately explained by stupidity.",
          "Robert Hanlon" },
        { "Developing the AI requires the work of a data scientist, and most of them understand "
          "neither data nor science.",
          "Scott LeGrand" },
        { "Before we work on artificial intelligence why don't we do something about natural "
          "stupidity?",
          "Steve Polyak" },
        { "If it weren't for C, we'd all be programming in BASI and OBOL.", "Anonymous" },
        { "Why is it that programmers always confuse Halloween with Christmas? Because 31 OCT = 25 "
          "DEC.",
          "Anonymous" },
        { "I'm not interrupting you, I'm putting our conversation in full-duplex mode.",
          "Antone Roundy" },
        { "The programmer got stuck in the shower because the instructions on the shampoo bottle "
          "said: Lather, Rinse, Repeat.",
          "Anonymous" },
        { "A programmer's spouse says 'While you're at the grocery store, buy some eggs.' The "
          "programmer never comes back.",
          "Anonymous" },
        { "What is a Unix or Linux sysadmin's favourite hangout place? Foo Bar.", "Anonymous" },
        { "If you have any trouble sounding condescending, find a UNIX user to show you how it's "
          "done.",
          "Scott Adams, Dilbert Cartoonist" },
        { "There is no place like ~", "Anonymous" },
        { "Thou shalt not kill -9", "Anonymous" },
        { "printf(\"%d is the year of the linux desktop\", year+1);", "Anonymous" },
        { "The use of COBOL cripples the mind; its teaching should therefore be regarded as a "
          "criminal offense.",
          "Edsger Dijkstra" },
        { "Computer system analysis is like child-rearing; you can do grievous damage, but you "
          "cannot ensure success.",
          "Tom DeMarcho" },
        { "PHP is a minor evil perpetrated and created by incompetent amateurs, whereas Perl is a "
          "great and insidious evil, perpetrated by skilled but perverted professionals.",
          "Jon Ribbens" },
        { "C is not a high-level language.", "Brian Kernighan, C author" },
        { "I will not be a lemming and follow the crowd over the cliff and into the C.",
          "John Beidler" },
        { "C is quirky, flawed, and an enormous success.", "Dennis Ritchie, C author" },
        { "A C program is like a fast dance on a newly waxed dance floor by people carrying "
          "razors.",
          "Waldi Ravens" },
        { "Fifty years of programming language research, and we end up with C++???",
          "Richard O'Keefe" },
        { "Quite frankly, even if the choice of C were to do *nothing* but keep the C++ "
          "programmers out, that in itself would be a huge reason to use C.",
          "Linus Torvalds" },
        { "Considering the current sad state of our computer programs, software development is "
          "clearly still a black art, and cannot yet be called an engineering discipline.",
          "William Jefferson Clinton" },
        { "I am rarely happier than when spending an entire day programming my computer to perform "
          "automatically a task that it would otherwise take me a good ten seconds to do by hand.",
          "Douglas Adams" },
        { "#define QUESTION ((bb),| !(bb))", "William Shakespeare" },
        { "I didn't know what MD was. I think I've managed to catch up.", "Berk Hess" },
        { "Teemu [Murtola] keeps beating our code, but that's fine because he's always right.",
          "Berk Hess" },
        { "Schrödinger's backup: The condition of any backup is unknown until a restore is "
          "attempted.",
          "Anonymous" },
        { "Don't waste pure thoughts on dirty enzymes.", "Efraim Racker" },
        { "I like single-molecule experiments because I hate to simulate 10^23 molecules at the "
          "same time.",
          "Helmut Grubmüller" },
        { "I wanted to make a clever chemistry joke, but the best ones Argon.", "39.948" },
        { "Not to get technical... but according to chemistry, alcohol is a solution.",
          "Anonymous" },
        { "The physical chemists never use their eyes and are most lamentably lacking in chemical "
          "culture. It is essential to cast out from our midst, root and branch, this physical "
          "element and return to our laboratories.",
          "Henry Edward Armstrong" },
        { "Time is the best appraiser of scientific work, and I am aware that an industrial "
          "discovery rarely produces all its fruit in the hands of its first inventor.",
          "Louis Pasteur" },
        { "Still I had a lurking question. Would it not be better if one could really 'see' "
          "whether molecules as complicated as the sterols, or strychnine were just as experiment "
          "suggested?",
          "Dorothy Hodgkin" },
        { "We think there is color, we think there is sweet, we think there is bitter, but in "
          "reality there are atoms and a void.",
          "Democritus" },
        { "A cop pulls Heisenberg over and asks him 'Do you know how fast you were going?' "
          "Heisenberg replies 'No, but I know exactly where I am'. The cop says 'You were doing 55 "
          "in a 35 zone'. Heisenberg: 'Great! Now I'm lost!",
          "Anonymous" },
        { "Two chemists walk into a bar. The first one says, 'I'll have some H2O.'. The second one "
          "says, 'I'll have some H2O, too'. He dies.",
          "Anonymous" },
        { "There are only two hard things in computer science - cache invalidation, naming things "
          "and off-by-one errors.",
          "Anonymous" },
        { "Science, for me, gives a partial explanation for life. In so far as it goes, it is "
          "based on fact, experience and experiment.",
          "Rosalind Franklin" },
        { "I was taught that the way of progress was neither swift nor easy.", "Marie Curie" },
        { "Life need not be easy, provided only that it is not empty.", "Lise Meitner" },
        { "We ignore public understanding of science at our peril.", "Eugenie Clark" },
        { "All sorts of things can happen when you're open to new ideas and playing around with "
          "things.",
          "Stephanie Kwolek, inventor of Kevlar" },
        { "As always in life, people want a simple answer... and it's always wrong.", "Marie Daly" },
        { "For a research worker the unforgotten moments of his life are those rare ones which "
          "come after years of plodding work, when the veil over natures secret seems suddenly to "
          "lift & when what was dark & chaotic appears in a clear & beautiful light & pattern.",
          "Gerty Cori" },
        { "The more clearly we can focus our attention on the wonders and realities of the "
          "universe about us, the less taste we shall have for destruction.",
          "Rachel Carson" },
        { "I didn't want to just know names of things. I remember really wanting to know how it "
          "all worked.",
          "Elizabeth Blackburn" },
        { "Science is not a boy's game, it's not a girl's game. It's everyone's game. It's about "
          "where we are and where we're going.",
          "Nichelle Nichols" },
        { "If you know you are on the right track, if you have this inner knowledge, then nobody "
          "can turn you off... no matter what they say.",
          "Barbara McClintock" },
        { "Science and everyday life cannot and should not be separated.", "Rosalind Franklin" },
        { "I hadn't been aware that there were doors closed to me until I started knocking on "
          "them.",
          "Gertrude Elion" },
        { "Humans are allergic to change. They love to say, 'We've always done it this way.' I try "
          "to fight that. That's why I have a clock on my wall that runs counter-clockwise.",
          "Grace Hopper" },
        { "Be less curious about people and more curious about ideas.", "Marie Curie" },
        { "I was a bit of an artist, and somewhere along the way had gotten the idea that "
          "computers could be used for animation and artists, because in-betweening was so "
          "tedious... Of course, everyone thought I was nuts.",
          "Carla Meninsky, Atari engineer" },
        { "I think it's very important to get more women into computing. My slogan is: Computing "
          "is too important to be left to men.",
          "Karen Jones" },
        { "There is only one thing worse than coming home from the lab to a sink full of dirty "
          "dishes, and that is not going to the lab at all!",
          "Chien-Shiung Wu" },
        { "They never asked me to go back over (my calculations) because when I did it, I had done "
          "my best, and it was right.",
          "Kathrine Jonson (legendary NASA mathematician)" },
        { "Science is a way of thinking much more than it is a body of knowledge.", "Carl Sagan" },
        { "Science is organized knowledge. Wisdom is organized life.", "Immanuel Kant" },
        { "Everything is theoretically impossible, until it is done.", "Robert Heinlein" },
        { "Bad times have a scientific value. These are occasions a good learner would not miss.",
          "Ralph Waldo Emerson" },
        { "The scientist is not a person who gives the right answers, he's one who asks the right "
          "questions.",
          "Claude Levi-Strauss" },
        { "Rockets are cool. There's no getting around that.", "Elon Musk" },
        { "Life would be tragic if it weren't funny.", "Stephen Hawking" },
        { "Somewhere, something incredible is waiting to be known.", "Carl Sagan" },
        { "I am driven by two main philosophies: know more today about the world than I knew "
          "yesterday and lessen the suffering of others. You'd be surprised how far that gets you.",
          "Neil deGrasse Tyson" },
        { "Your assumptions are your windows on the world. Scrub them off every once in a while, "
          "or the light won't come in.",
          "Isaac Asimov" },
        { "Religion is a culture of faith; science is a culture of doubt.", "Richard Feynman" },
        { "Millions saw the apple fall, Newton was the only one who asked why?", "Bernard Baruch" },
        { "If you thought that science was certain - well, that is just an error on your part.",
          "Richard Feynman" },
        { "When my information changes, I alter my conclusions. What do you do, sir?",
          "John Maynard Keynes" },
        { "The aim of science is not to open the door to infinite wisdom, but to set a limit to "
          "infinite error.",
          "Bertolt Brecht, Life of Galileo" },
        { "As we all know, blinking lights means science.", "Joss Whedon" },
        { "But in my opinion, all things in nature occur mathematically.", "Rene Decartes" },
        { "A weed scientist goes into a shop. He asks: 'Hey, you got any of that inhibitor of "
          "3-phosphoshikimate-carboxyvinyl transferase?' Shopkeeper: 'You mean Roundup?' "
          "Scientist: 'Yeah, that's it. I can never remember that dang name!'",
          "John Pickett" },
        { "It is not clear that intelligence has any long-term survival value.", "Stephen Hawking" },
        { "The greatest shortcoming of the human race is our inability to understand the "
          "exponential function.",
          "Albert Bartlett" },
        { "You can get into a habit of thought in which you enjoy making fun of all those other "
          "people who don't see things as clearly as you do. We have to guard carefully against "
          "it.",
          "Carl Sagan" },
        { "I have no responsibility to live up to what others expect of me. That's their mistake, "
          "not my failing.",
          "Richard Feynman" },
        { "Highly organized research is guaranteed to produce nothing new.", "Frank Herbert" },
        { "Those who cannot remember the past are condemned to compute it.", "Steve Pinker" },
        { "If a rat is a good model for your emotional life, you're in big trouble.",
          "Robert Sapolsky" },
        { "I don't know how many of you have ever met Dijkstra, but you probably know that "
          "arrogance in computer science is measured in nano-Dijkstras.",
          "Alan Kay" },
        { "NASA spent millions of dollars inventing the ball-point pen so they could write in "
          "space. The Russians took a pencil.",
          "Will Chabot" },
        { "By denying scientific principles, one may maintain any paradox.", "Galileo Galilei" },
        { "Perfect is the enemy of good.", "Voltaire" },
        { "Men love to wonder, and that is the seed of science.", "Ralph Waldo Emerson" },
        { "In mathematics you don't understand things, you just get used to them",
          "John von Neumann" },
        { "A real scientist solves problems, not wails that they are unsolvable.",
          "Anne McCaffrey" },
        { "Science progresses best when observations force us to alter our preconceptions.",
          "Vera Rubin" },
        { "Our two greatest problems are gravity and paper work. We can lick gravity, but "
          "sometimes the paperwork is overwhelming.",
          "Wernher von Braun" },
        { "I have had my results for a long time, but I do not yet know how I am to arrive at "
          "them.",
          "Carl Friedrich Gauss" },
        { "The difficulty lies, not in the new ideas, but in escaping the old ones.",
          "John Maynard Keynes" },
        { "The way to succeed is to double your failure rate.", "Thomas J. Watson" },
        { "We are continually faced by great opportunities brilliantly disguised as insoluble "
          "problems.",
          "Lee Iacocca" },
        { "Mathematics is a game played according to certain rules with meaningless marks on "
          "paper.",
          "David Hilbert" },
        { "Mathematics is no more computation than typing is literature.", "John Allen Paulos" },
        { "Mathematics is like love; a simple idea, but it can get complicated.", "Anonymous" },
        { "I had a polynomial once. My doctor removed it.", "Michael Grant" },
        { "Since the mathematicians have invaded the theory of relativity I do not understand it "
          "myself any more.",
          "Albert Einstein" },
        { "I couldn't claim that I was smarter than sixty-five other guys - but the average of "
          "sixty-five other guys, certainly!",
          "Richard Feynman" },
        { "Your Excellency, I have no need of this hypothesis.",
          "Pierre Laplace, to Napoleon on why his works on celestial mechanics make no mention of "
          "God." },
        { "The sign of wisdom is to have more questions than answers.", "Abhijit Naskar" },
        { "Stupidity got us into this mess, and stupidity will get us out.", "Homer Simpson" },
        { "Trying is the first step towards failure.", "Homer Simpson" },
        { "Same sex marriage is not a gay privilege, it's equal rights. Privilege would be "
          "something like gay people not paying taxes. Like churches don't.",
          "Ricky Gervais" },
        { "Remember, being healthy is basically dying as slowly as possible.", "Ricky Gervais" },
        { "Pain is inevitable. Suffering is optional.", "Haruki Murakami" },
        { "Aristotle maintained that women have fewer teeth than men; although he was twice "
          "married, it never occurred to him to verify this statement by examining his wives' "
          "mouths.",
          "Bertrand Russell" },
        { "I don't believe in astrology; I'm a Sagittarian and we're skeptical.",
          "Arthur C. Clarke" },
        { "I see they found out the universe is 80 million years older than we thought. It's also "
          "been lying about its weight.",
          "Bill Maher" },
        { "Do you have mole problems? If so, call Avogadro at 602-1023.", "Jay Leno" },
        { "Marie, you're looking more radiant every day!", "Pierre Curie" },
        { "I don't want to achieve immortality through my work... I want to achieve it through not "
          "dying!",
          "Woody Allen" },
        { "Well, I am a dilettante. It's only in England that dilettantism is considered a bad "
          "thing. In other countries it's called interdisciplinary research.",
          "Brian Eno" },
        { "I think it would be a good idea.",
          "Mahatma Gandhi, when asked what he thought of Western civilization" },
        { "Nobody ever complained a seminar was too easy to understand.", "Ken Dill" },
        { "Academia is kind of like applied Marxism. The workers really do own the means of "
          "production.",
          "Niklas Blomberg" },
        { "The Lord of the Rings can be confusing to follow because many of the bad minions look "
          "and sound familiar; that's why Tolkien gave them each an ORCid.",
          "Caroline Bartman" },
        { "Mendeleev's first attempt, the perfluoric table, was a total disaster, and his "
          "subsequent attempts, the perchloric and perbromic tables, were not favorably received. "
          "Only his fourth attempt, the periodic table, gained general acceptance.",
          "Anonymous" },
        { "Don’t bring an anecdote to a data fight.", "Molly Hodgdon" },
        { "Give someone a program, you frustrate them for a day; teach them how to program, you "
          "frustrate them for a lifetime.",
          "David Leinweber" },
        { "Or (horrors!) use Berendsen!", "Justin Lemkul" },
        { "The absence of real intelligence doesn't prove you're using AI", "Magnus Lundborg" },
        { "People who do QM/MM must be rather patient and enjoy quality over speed",
          "Kresten Lindorff-Larsen" },
        { "I don’t think we’re afraid of inline assembly.", "Szilard Pall" },
        { "I'm a strong believer that ignorance is important in science. If you know too much, you "
          "start seeing reasons why things won't work. That's why its important to change your "
          "field to collect more ignorance.",
          "Sydney Brenner" },
        { "It's more useful when you know what you're doing.", "Artem Zhmurov" },
        { "I have noticed a large, negative correlation between having a well-defined mission "
          "workload and concern for the Top500. It's almost like LINPACK is what you focus on when "
          "you don't know what to focus on.",
          "Jeff Hammond" },
        { "Between equal rights, force decides.", "Karl Marx" },
        { "To dissimulate is to feign not to have what one has. To simulate is to feign to have "
          "what one hasn't.",
          "Jean Baudrillard" },
        { "Install our Free Energy Patents app! There is energy all around us; and it's free! "
          "Free energy is everywhere, and all around you, just waiting to be extracted! Over "
          "100+ free energy patents!",
          "Mind and Miracle Productions on Twitter, spamming a FEP thread" },
        { "\"A slow sort of country!\" said the Queen. \"Now, HERE, you see, it "
          "takes all the running YOU can do, to keep in the same place. If you want "
          "to get somewhere else, you must run at least twice as fast as that!\"",
          "Lewis Carroll" },
        { "More than 10000000 total errors detected.  I'm not reporting any more. "
          "Final error counts will be inaccurate.  Go fix your program!",
          "Valgrind while memory debugging mdrun" },
        { "If we are going to have SYCL, can we have a hammer as well?", "Joe Jordan" },
        { "We can make it into a friend class. But I don't like having friends.", "Joe Jordan" },
        { "A method is more important than a discovery, since the right method will lead to new "
          "and even more important discoveries.",
          "Lev Landau" },
        { "Product of optimism and knowledge is a constant.", "Lev Landau" },
        { "Why add prime numbers? Prime numbers are made to be multiplied.", "Lev Landau" },
        { "How wonderful that we have met with a paradox. Now we have some hope of making "
          "progress.",
          "Niels Bohr" },
        { "We must be clear that when it comes to atoms, language can be used only as in poetry. ",
          "Niels Bohr" },
        { "\"What are the biological implications of your research?\" - \"Well, I simulate "
          "water.\" ",
          "Petter Johansson" },
        { "Everything what mathematicians were saying for the last 50 years is slowly catching up "
          "with us.",
          "David van der Spoel" },
        { "I tend to consider myself as a scientist.",
          "Emmanuelle Charpentier, when asked about the importance of two women sharing the Nobel "
          "Prize for Chemistry" },
        { "I identified myself very early on as a scientist rather than a student - as someone "
          "creating knowledge rather than simply absorbing it.",
          "Emmanuelle Charpentier" },
        { "Look, I don't want to compete, so let's divide up physics between us. I'll take auroras "
          "and you take the rest of the universe.",
          "Joan Feynman to her brother Richard" },
        { "There are three kinds of men. The one that learns by reading. The few who learn by "
          "observation. The rest of them have to pee on the electric fence for themselves.",
          "Will Rogers" },
        { "I can't help but think the model is ungrateful for all that nice data I gave it. Jerk.",
          "Kate Stafford" },
        { "What do you call an acid with an attitude? A-mean-oh acid.", "Anonymous" },
        { "Science grows like a weed every year.", "Kary Mullis" },
        { "A good sign these days when you're listening to someone talk about this epidemic is the "
          "number of times they say 'We don't know yet'. The more of those, the better.",
          "Derek Lowe" },
        { "Nullis in verba [Nobody's word is final].", "Motto of the Royal Society" },
        { "Calling a system 'non-linear' is like calling all wild animals 'non-elephants'.",
          "Stan Ulam" },
        { "Given enough eyeballs, all bugs are shallow.",
          "Linus Torvalds, on the power of open source" },
        { "If every study was groundbreaking, we'd end up with a bunch of holes in the ground and "
          "nothing built.",
          "Anonymous" },
        { "Expertise is not inherently good.", "Joe Jordan" },
        { "I couldn't give a shit about ribosomes.",
          "Björn Forsberg, presenting his thesis, including two papers on ribosomes" },
        { "Here's something fun to do: Next time you approach a conversation say 'I want to talk "
          "to someone technical... Oh! There's a woman!' and walk straight over to her.",
          "Patty Lopez" },
        { "After a few talks we usually sit down to do some work... or drinking.", "Mike Klein" },
        { "The message here is that Thermodynamic Integration sucks.", "Berk Hess" },
        { "Enthusiasm is the mother of effort, and without it nothing great was ever achieved.",
          "Ralph Waldo Emerson" },
        { "Our hands are tied by physics.", "Christian Blau" },
        { "If all it takes to motivate you is a fancy picture and quote, you probably have a very "
          "easy job. The type of job computers will soon be doing.",
          "Anonymous" },
        { "At school I had a teacher that didn't like me and I didn't like him. At the end of the "
          "year he decided to fail me. The ironic thing is that the topic was chemistry. I have "
          "the distinction of being the only Chemistry Laurate who failed the topic in high "
          "school.",
          "Thomas Lindahl" },
        { "Chemistry: It tends to be a messy science.",
          "Gunnar von Heijne, former chair of the Nobel Committee for chemistry" },
        { "Computers are incredibly fast, accurate and stupid. Humans are incredibly slow, "
          "inaccurate and... also stupid.",
          "Anonymous" },
        { "Schrödinger's backup: The condition of any backup is unknown until a restore is "
          "attempted.",
          "Anonymous" },
        { "If my PhD doesn't allow me to be right on the internet, what is it even good for?",
          "Martin Vögele" },
        { "A little less conversation, a little more action, please.", "Elvis Presley" },
        { "Friends don't let friends use Berendsen!", "John Chodera (on Twitter)" },
        { "The plural of regex is regrets", "Steve McCarthy (on Twitter)" },
        { "Culture eats strategy for breakfast", "Peter Drucker" },
        { "Roses are read // Violets are blue // Unexpected '}' on line 32", "Anonymous" },
        { "We cannot wait for Nature's good graces - to take them from her is our goal",
          "Ivan Michurin" },
        { "By three methods we may learn wisdom: First, by reflection, which is noblest; "
          "Second, by imitation, which is easiest; "
          "and third by experience, which is the bitterest.",
          "Confucius" },
        { "There are three types of people: Those who see, those who see when they are shown, "
          "and those who do not see.",
          "Leonardo da Vinci" },
        { "Hey, it's me - Pandora. Welcome to my new unboxing video!", "Anonymous" },
        { "It is an unfortunate fact that when you raise the question of the reliability of many "
          "simulations you are often told about how much manpower went into it, how large & fast "
          "the computer is, how important the problem is, and such things, which are completely "
          "irrelevant to the question that was asked.",
          "Richard Hamming" },
        { "Sie haben also recht gehabt, Sie Spitzbube. [You were right after all, you rascal.]",
          "Albert Einstein (letter to Wolfgang Pauli)" },
        { "You should call it 'entropy'. No one knows what entropy really is, so in a debate you "
          "will always have the advantage.",
          "John von Neumann to Claude Shannon, "
          "on why he should borrow the term for information theory" },
        { "We must have perseverance and above all confidence in ourselves. We must believe that "
          "we are gifted for something and that this thing must be attained.",
          "Marie Curie" },
        { "Courage is like - it's a habitus, a habit, a virtue: you get it by courageous acts. "
          "It's like you learn to swim by swimming. You learn courage by couraging.",
          "Marie Daly" },
        { "We will always have STEM with us. Some things will drop out of the public eye and "
          "will go away, but there will always be science, engineering, and technology. And there "
          "will always, always be mathematics.",
          "Katherine Johnson" },
        { "If you want to change the future, start living as if you're already there.",
          "Lynn Conway" },
        { "What you do makes a difference, and you have to decide what kind of difference you want "
          "to make.",
          "Jane Goodall" },
        { "I don't fear death because I don't fear anything I don't understand.", "Hedy Lamarr" },
        { "You cannot hope to build a better world without improving the individuals.",
          "Marie Curie" },
        { "Forget this world and all its troubles and, if possible, it's multitudinous Charlatans "
          "- everything, in short, but the Enchantress of Numbers",
          "Ada Lovelace" },
        { "We had the quaint notion at the time that software should be completely, absolutely "
          "free of bugs. Unfortunately it's a notion that never really quite caught on.",
          "Mary Allen Wilkes" },
        { "Like what you do, and then you will do your best.", "Katherine Johnson" },
        { "All creative people want to do the unexpected", "Hedy Lamarr" },
        { "Unlike teachers or doctors, our efforts improve the lives of people we'll never meet.",
          "Katie Busch-Sorensen" },
        { "Numbers have life; they’re not just symbols on paper.", "Shakuntala Devi" },
        { "Not everyone is capable of madness; and of those lucky enough to be capable, not many "
          "have the courage for it.",
          "August Strindberg" },
        { "The historical meaning of QM/MM: Quanto Mangi, Mamma Mia!",
          "Attilio Vittorio Vargiu, at BioExcel 2022 Summer School in Sardinia" },
        { "I wanted to have a real language where programmers could write real programs, "
          "and see whether they would find the idea of data abstraction at all useful.",
          "Barbara Liskov" },
        { "I never thought about breaking barriers, I was just interested in what I was doing "
          "and kept going.",
          "Barbara Liskov" },
        { "My husband was on the Internet everyday when I got the Turing award, and one day "
          "he saw a quote by someone who said 'Why did she get the award? Everyone already knows "
          "this!'",
          "Barbara Liskov" },
        { "When I got my first job as a programmer there were quite a number of women doing "
          "programming "
          "because there was no computer science education and people were hired from many "
          "different "
          "fields if it seemed like they could do the work. It was better then, probably, in terms "
          "of "
          "proportions, not necessarily in women being paid as much as men and so forth, but "
          "definitely, "
          "there were women around.",
          "Barbara Liskov" },
        { "Prior to 1965 there were none, and after 1965 there was a nun.",
          "Sister Mary Kenneth Keller regarding women with PhDs in computer science" },
        { "We're having an information explosion, among others, and it's certainly obvious that "
          "information is of no use unless it's available.",
          "Sister Mary Kenneth Keller" },
        { "For the first time we can now mechanically simulate the cognitive process. We can "
          "make studies in artificial intelligence. Beyond that, this mechanism can be used to "
          "assist "
          "humans in learning. As we are going to have more mature students in greater numbers as "
          "time "
          "goes on, this type of teaching will probably be increasingly important.",
          "Sister Mary Kenneth Keller" },
        { "Why would I pay GitHub $100/year for an AI to tell me what code to write when men "
          "do it for free?",
          "Safia Abdalla" },
        { "So can we think of the canonical ensemble sort of like if you ran the same "
          "season of temptation island over and over again and took statistics about how many "
          "times each person hooks up with each other person?",
          "Teddy Press" },
        { "The biggest lie in science is 'data and code available upon request'", "Michael Eisen" },
        { "What do you want out of life?", "Jack Kerouac, On The Road" },
        { "Developer accused of unreadable code refuses to comment", "Molly Struve" },
        { "It's hard to ignore 12 orders of magnitude", "John Urbanic" },
        { "Success is going from one failure to another without loss of enthusiasm",
          "Winston Churchill" },
        { "Apologies to the astrophysics student I met at a party years ago. When you told me "
          "how many hours a day you used 4chan and how much you love it, I gave you a funny "
          "look and walked away. Now, a decade later, I realize you were talking about Fortran.",
          "Anonymous" },
        { "The only place success comes before work is in the dictionary", "Vince Lombardi" },
        { "They call me 007: 0 pull requests reviewed, 0 features completed, 7 bugs created.",
          "Anonymous" },
        { "Optimist: The glass is 1/2 full. Pessimist: The glass is 1/2 empty. "
          "Excel: The glass is January 2nd.",
          "John Feminella" },
        { "FORTRAN. Input: reason, output: pleasure", "ORDA, FORTRAN board game" },
        { "gmx fellowship-writing -g grant_name -s protein_structure_involved -o output -m "
          "method_used -p list_of_pi",
          "Tanadet Pipatpolkai, while discussing new features for GROMACS" },
        { "I came up with the new convergence method, it's called a deadline driven convergence. "
          "My simulation is converged when it hits the deadline.",
          "Tanadet Pipatpolkai" },
        { "Lets get back to beer", "Yuxuan Zhuang, in a discussion about science communication" },
        { "You ONLY have to do the coding ...",
          "Anton Jansen, to core developer, on implementing new features" },
        { "There are way too many quotes", "Sebastian Wingbermuehle" },
        { "It is not critical to add the next quote to a patch release", "Paul Bauer" },
        { "It is a cute toxin.", "Rebecca Howard" },
        { "Everything is failing", "Paul Bauer" },
        { "Requiem, bring the dissident from slumber", "Bad Religion" },
        { "I can't relate to you", "Bad Religion" },
        { "You are wrong!", "NOFX" },
        { "The final page is written in the books of history", "Bad Religion" },
        { "Would you give it all up to live again?", "Bad Religion" },

    };

    if (beCool())
    {
        auto quote = getPseudoRandomElement<Quote>(quoteArray);
        return formatString("GROMACS reminds you: \"%s\" (%s)", quote.text, quote.author);
    }
    else
    {
        return "Thanx for Using GROMACS - Have a Nice Day";
    }
}

} // namespace gmx
