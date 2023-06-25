package Utilities;

public class Table {

    private int columns;
    public Object[][] items;

    /**
     * Constructs a table with the specified number of columns.
     *
     * @param columns the number of columns in the table
     */
    public Table(int columns) {
        this.columns = columns;
        items = new Object[0][columns];
    }

    /**
     * Adds a new row to the table.
     *
     * @param row the row to be added
     */
    public void addNewRow(Object[] row) {
        Object[][] tempItems = items;
        items = new Object[items.length + 1][columns];
        System.arraycopy(tempItems, 0, items, 0, tempItems.length);
        items[items.length - 1] = row;
    }

    /**
     * Prints all the items in the table.
     */
    public void printItems() {
        for (Object[] objs : items) {
            for (Object obj : objs) {
                System.out.print(obj + " ; ");
            }
            System.out.println();
        }
    }

    /**
     * Edits a specific cell in the table.
     *
     * @param rowIndex the index of the row to be edited
     * @param columnIndex the index of the column to be edited
     * @param newData the new data to be set in the cell
     */
    public void editCell(int rowIndex, int columnIndex, Object newData) {
        items[rowIndex][columnIndex] = newData;
    }

    /**
     * Merges two arrays into a single array.
     *
     * @param array1 the first array to be merged
     * @param array2 the second array to be merged
     * @return the merged array
     */
    public Object[] mergeArrays(Object[] array1, Object[] array2) {
        int a1 = array1.length;
        int a2 = array2.length;
        Object[] mergedArray = new Object[a1 + a2];
        System.arraycopy(array1, 0, mergedArray, 0, a1);
        System.arraycopy(array2, 0, mergedArray, a1, a2);
        return mergedArray;
    }

}
