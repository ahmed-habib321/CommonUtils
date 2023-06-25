package Utilities;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;
import javax.imageio.ImageIO;
import javax.swing.*;

public class Tools {

    /**
     * Displays a message dialog box with the given message.
     *
     * @param message The message to display.
     */
    public static void showMessageDialog(String message) {
        JOptionPane.showMessageDialog(null, message);
    }

    /**
     * Displays a confirmation dialog box with the given message.
     *
     * @param message The message to display.
     * @return true if the user selects "Yes," false otherwise.
     */
    public static boolean showConfirmationDialog(String message) {
        int result = JOptionPane.showConfirmDialog(null, message);
        return result == JOptionPane.YES_OPTION;
    }

    /**
     * Creates a new folder with the given name at the specified path.
     *
     * @param folderName The name of the folder to create.
     * @param path The path where the folder should be created.
     */
    public static void createFolder(String folderName, String path) {
        File folder = new File(path + "/" + folderName);
        folder.mkdir();
    }

    /**
     * Creates a new folder with the given name.
     *
     * @param folderName The name of the folder to create.
     */
    public static void createFolder(String folderName) {
        File folder = new File(folderName);
        folder.mkdir();
    }

    /**
     * Opens a JFrame form, sets its properties, and displays it.
     *
     * @param form The JFrame form to open.
     */
    public static void openForm(JFrame form) {
        try {
            form.setLocationRelativeTo(null);
            Image icon = ImageIO.read(Tools.class.getResource("emp.png"));
            form.setIconImage(icon);
            form.setTitle("Company");
            form.getContentPane().setBackground(Color.WHITE);
            form.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            form.setVisible(true);
            form.setResizable(false);
        } catch (IOException ex) {
            showMessageDialog(ex.getMessage());
        }
    }

    /**
     * Clears the text of all JTextFields in the given container and its nested
     * containers.
     *
     * @param container The container to clear.
     */
    public static void clearTextFields(Container container) {
        for (Component component : container.getComponents()) {
            if (component instanceof JTextField) {
                JTextField textField = (JTextField) component;
                textField.setText("");
            } else if (component instanceof Container) {
                clearTextFields((Container) component);
            }
        }
    }

    /**
     * Creates an empty file with the given file name.
     *
     * @param fileName The name of the file to create.
     */
    public static void createEmptyFile(String fileName) {
        try {
            File file = new File(fileName);
            file.createNewFile();
        } catch (IOException ex) {
            showMessageDialog(ex.getMessage());
        }
    }

    /**
     * Creates empty files with the given file names.
     *
     * @param fileNames The names of the files to create.
     */
    public static void createEmptyFiles(String[] fileNames) {
        for (String fileName : fileNames) {
            createEmptyFile(fileName);
        }
    }

    /**
     * Creates a file with the given file name and writes the data to it.
     *
     * @param fileName The name of the file to create.
     * @param data The data to write to the file.
     */
    public static void createFile(String fileName, Object[] data) {
        try {
            PrintWriter writer = new PrintWriter(fileName);
            for (Object obj : data) {
                writer.println(obj);
            }
            writer.close();
        } catch (FileNotFoundException ex) {
            showMessageDialog(ex.getMessage());
        }
    }

    /**
     * Creates files with the given file names and writes the corresponding data
     * to each file.
     *
     * @param fileNames The names of the files to create.
     * @param allData The data to write to the files.
     */
    public static void createFiles(String[] fileNames, Object[][] allData) {
        for (int i = 0; i < fileNames.length; i++) {
            createFile(fileNames[i], allData[i]);
        }
    }

    /**
     * Displays an input dialog box and returns the user's input as an Object.
     *
     * @param title The title of the input dialog box.
     * @return The user's input as an Object.
     */
    public static Object showInputDialog(String title) {
        return JOptionPane.showInputDialog(title);
    }

    /**
     * Extracts and returns only the digits from the given text.
     *
     * @param text The text to extract digits from.
     * @return A string containing only the digits.
     */
    public static String extractDigits(String text) {
        StringBuilder digits = new StringBuilder();
        for (char c : text.toCharArray()) {
            if (Character.isDigit(c)) {
                digits.append(c);
            }
        }
        return digits.toString();
    }

    /**
     * Extracts and returns the integer value represented by the digits in the
     * given text.
     *
     * @param text The text to extract the integer from.
     * @return The integer value represented by the digits.
     */
    public static int extractInteger(String text) {
        String digits = extractDigits(text);
        return Integer.parseInt(digits);
    }

    /**
     * Removes all digits from the given text and returns the result.
     *
     * @param text The text to remove digits from.
     * @return The text without any digits.
     */
    public static String removeDigits(String text) {
        StringBuilder result = new StringBuilder();
        for (char c : text.toCharArray()) {
            if (!Character.isDigit(c)) {
                result.append(c);
            }
        }
        return result.toString();
    }

    /**
     * Takes a screenshot of the entire screen and saves it as an image with the
     * given name.
     *
     * @param imageName The name of the image file to save.
     */
    public static void captureScreen(String imageName) {
        try {
            Robot robot = new Robot();
            Rectangle screenRect = new Rectangle(Toolkit.getDefaultToolkit().getScreenSize());
            BufferedImage screenshot = robot.createScreenCapture(screenRect);
            ImageIO.write(screenshot, "jpg", new File(imageName + ".jpg"));
        } catch (AWTException | HeadlessException | IOException ex) {
            showMessageDialog(ex.getMessage());
        }
    }

    /**
     * Takes a screenshot of the specified JFrame form, saves it as an image
     * with the given name, and restores the form to its previous state.
     *
     * @param imageName The name of the image file to save.
     * @param form The JFrame form to take a screenshot of.
     */
    public static void captureScreen(String imageName, JFrame form) {
        try {
            form.setState(JFrame.ICONIFIED);
            Robot robot = new Robot();
            Rectangle screenRect = new Rectangle(Toolkit.getDefaultToolkit().getScreenSize());
            BufferedImage screenshot = robot.createScreenCapture(screenRect);
            ImageIO.write(screenshot, "jpg", new File(imageName + ".jpg"));
            form.setState(JFrame.NORMAL);
        } catch (AWTException | HeadlessException | IOException ex) {
            showMessageDialog(ex.getMessage());
        }
    }

    /**
     * Sets the appearance of a JTable's header.
     *
     * @param table The JTable to modify.
     */
    public static void setTableHeaderAppearance(JTable table) {
        table.getTableHeader().setBackground(new Color(0, 0, 100));
        table.getTableHeader().setForeground(Color.WHITE);
    }
}
