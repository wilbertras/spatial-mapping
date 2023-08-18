# Python 3 program of the above approach
 
# Stores the 8 possible combinations of
# moves that the knight can follow
DirX = [2, 1, -1, -2, -2, -1, 1, 2]
DirY = [1, 2, 2, 1, -1, -2, -2, -1]
 
# Function to find if (i, j) is a valid
# cell for the knight to move and it
# exists within the chessboard
def isSafe(i, j, n, Board):
    return i >= 0 and j >= 0 and i < n and j < n and Board[i][j] == 0
 
# Stores whether there exist any valid path
isPossible = False
# Recursive function to iterate through all
# the paths that the knight can follow
def knightTour(ChessBoard, N, x, y, visited=1):
    global isPossible
    n = 0
    req = 1
    # Mark the current square of the chessboard
    ChessBoard[x][y] = visited
 
    # If the number of visited squares are equal
    # to the total number of squares
    if visited == N * N:
        isPossible = True
 
        # Print the current state of ChessBoard
        print(ChessBoard)
        
        n += 1
        # # Backtrack to the previous move
        # ChessBoard[x][y] = 0
        # return

    if n < req:
        # Iterate through all the eight possible moves
        # for a knight
        for i in range(8):

            # Stores the new position of the knight
            # after a move
            newX = x + DirX[i]
            newY = y + DirY[i]

            # If the new position is a valid position
            # recursively call for the next move
            if isSafe(newX, newY, N, ChessBoard) and not ChessBoard[newX][newY]:
                knightTour(ChessBoard, N, newX, newY, visited + 1)
        
        # Backtrack to the previous move
        ChessBoard[x][y] = 0


def init(N, X, Y):
    ChessBoard = [[0 for j in range(N)] for i in range(N)]
    N = len(ChessBoard)
    
    knightTour(ChessBoard, N, X - 1, Y - 1)

    # If no valid sequence of moves exist
    if not isPossible:
        print(-1)