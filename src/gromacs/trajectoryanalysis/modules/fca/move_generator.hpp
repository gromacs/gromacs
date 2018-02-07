#ifndef MOVE_GENERATOR_HPP
#define MOVE_GENERATOR_HPP

#include <vector>
#include "utils.hpp"

class MoveGenerator{
    static real defaultMin(){
        return real(-1e15);
    }

    struct Ijval{
        Ijval(): other(-1), val(defaultMin()), inOrder(false){
        }
        int other;
        real val;
        bool inOrder;
    };

    const int dim;
    std::vector<Ijval> moves;
    int currentMinMoveIndex;

    const int kick_max;
    // if moves get replaced due to axis overlap, they might be possible later on --- store them in kickouts
    std::vector<std::pair<int,int>> kickouts;

public:
    explicit MoveGenerator(const int inDim, const int inMaxKickValues)
        : dim(inDim),
          currentMinMoveIndex(-1),
          kick_max(inMaxKickValues){
        moves.resize(dim);
        kickouts.reserve(kick_max);
    }

    real getActualSum() const{
        real sum_val      = 0;
        for(const Ijval& ijval : moves){
            sum_val += ijval.val;
        }
        return sum_val/2;
    }

    int getNumberOfKicks() const{
        return static_cast<int>(kickouts.size());
    }

    int getActualNumberOfMoves() const{
        int nbMoves = 0;
        for(const Ijval& ijval : moves){
            if(ijval.other != -1){
                nbMoves += 1;
            }
        }
        return nbMoves;
    }

    std::vector<std::pair<int,int>> getMoves() const {
        std::vector<std::pair<int,int>> movesPairs;
        movesPairs.reserve(moves.size());
        for(int idx = 0 ; idx < static_cast<int>(moves.size()) ; ++idx){
            const Ijval& ijval = moves[idx];
            if(idx < ijval.other){
                if(ijval.inOrder){
                    movesPairs.emplace_back(idx, ijval.other);
                }
                else{
                    movesPairs.emplace_back(ijval.other, idx);
                }
            }
        }
        return movesPairs;
    }

    real getCurrentMin() const{
        if(currentMinMoveIndex == -1){
            return defaultMin();
        }
        else{
            return moves[currentMinMoveIndex].val;
        }
    }

    int getCurrentMinIdx() const{
        return currentMinMoveIndex;
    }

    std::vector<std::pair<int,int>> getKickPairs() const{
        return kickouts;
    }

    void add(const int i, const int j, const real dum){
        if(moves[i].other == -1 && moves[j].other == -1){
            moves[i].other = j;
            moves[i].val = dum;
            moves[i].inOrder = (i <= j);
            moves[j].other = i;
            moves[j].val = dum;
            moves[j].inOrder = (i <= j);
            if(currentMinMoveIndex == -1 || moves[currentMinMoveIndex].val > moves[i].val){
                currentMinMoveIndex = i;
            }
        }
        else {
            real testval = 0;
            if(moves[i].other != -1){
                testval += moves[i].val;
            }
            if(moves[j].other != -1){
                testval += moves[j].val;
            }
            if(testval < dum){
                if(moves[i].other != -1){
                    if(static_cast<int>(kickouts.size()) < kick_max){
                        kickouts.emplace_back(i, moves[i].other);
                    }
                    assert(moves[moves[i].other].other == i);
                    moves[moves[i].other].other = -1;
                    moves[moves[i].other].val = defaultMin();
                    currentMinMoveIndex = moves[i].other;
                }
                if(moves[j].other != -1 && moves[j].other != i){
                    if(static_cast<int>(kickouts.size()) < kick_max){
                        kickouts.emplace_back(j, moves[j].other);
                    }
                    assert(moves[moves[j].other].other == j);
                    moves[moves[j].other].other = -1;
                    moves[moves[j].other].val = defaultMin();
                    currentMinMoveIndex = moves[j].other;
                }
                moves[i].other = j;
                moves[i].val = dum;
                moves[i].inOrder = (i <= j);
                moves[j].other = i;
                moves[j].val = dum;
                moves[j].inOrder = (i <= j);
                if(currentMinMoveIndex == -1 || moves[currentMinMoveIndex].val > moves[i].val){
                    currentMinMoveIndex = i;
                }
            }
            else if(static_cast<int>(kickouts.size()) < kick_max){
                kickouts.emplace_back(i, j);
            }
        }
    }
};


#endif
