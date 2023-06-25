package Utilities;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;
import javax.swing.JTable;
import javax.swing.table.DefaultTableModel;

/**
 * A utility class for saving and loading data from a JTable with two columns.
 */
public class SaveLoadUtility {

    /**
     * Saves the data from the JTable into a file.
     *
     * @param table The JTable containing the data to be saved.
     * @param filename The name of the file to save the data into.
     */
    public static void save(JTable table, String filename) {
        DefaultTableModel model = (DefaultTableModel) table.getModel();
        Vector<Vector> dataVector = model.getDataVector();

        try ( PrintWriter saveFile = new PrintWriter(filename + ".txt")) {
            for (Vector vector : dataVector) {
                String save = vector.get(0) + ";" + vector.get(1);
                saveFile.println(save);
            }
        } catch (IOException ex) {
            // Handle the exception here or log it for debugging.
            ex.printStackTrace();
        }
    }

    /**
     * Loads data from a file into the JTable.
     *
     * @param table The JTable to load the data into.
     * @param filename The name of the file to load the data from.
     */
    public static void load(JTable table, String filename) {
        DefaultTableModel model = (DefaultTableModel) table.getModel();

        String filePath = filename + ".txt";

        try ( BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                model.addRow(line.split(";"));
            }
        } catch (IOException e) {
            // Handle the exception here or log it for debugging.
            e.printStackTrace();
        }
    }
}
