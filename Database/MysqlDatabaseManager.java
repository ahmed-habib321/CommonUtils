package Database;

import java.sql.*;
import java.util.Vector;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JTable;
import javax.swing.table.DefaultTableModel;
import CustomGUI.CustomComboBox;
import Utilities.Tools;
import Utilities.Table;

/**
 * This class provides methods for interacting with a database using JDBC.
 */
public class MysqlDatabaseManager {

    private static String URL = "";
    private static final String USER = "root";
    private static final String PASS = "";
    private static String DBName = "";

    private static Connection con;

    /**
     * Sets the database URL.
     */
    private static void setURL() {
        URL = "jdbc:mysql://localhost:3306/" + DBName;

    }

    /**
     * Sets up the database connection.
     */
    private static void setConnection() {
        try {
            //Class.forName("oracle.jdbc.driver.OracleDriver");
            setURL();
            con = DriverManager.getConnection(URL, USER, PASS);
        } catch (SQLException ex) {
            Tools.showMessageDialog(ex.getMessage());
        }
    }

    /**
     * Checks if the provided username and password are valid for login.
     *
     * @param username the username to check
     * @param password the password to check
     * @return true if the login is successful, false otherwise
     */
    public static boolean checkLogin(String username, String password) {
        try {
            setConnection();
            Statement st = con.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
            String check = "SELECT * FROM users WHERE username='" + username + "' AND password='" + password + "'";
            st.executeQuery(check);
            while (st.getResultSet().next()) {
                con.close();
                //tools.tools.msgBox("ok");
                return true;
            }
            con.close();
        } catch (SQLException e) {
            Tools.showMessageDialog(e.getMessage());
        }
        return false;
    }

    /**
     * Executes an SQL statement for editing data in the database.
     *
     * @param statement the SQL statement to execute
     * @return true if the execution is successful, false otherwise
     */
    public static boolean doEdit(String statement) {
        try {
            setConnection();
            Statement st = con.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
            st.execute(statement);
            con.close();
            return true;
        } catch (SQLException ex) {
            Tools.showMessageDialog(ex.getMessage());
            try {
                con.close();
            } catch (SQLException ex1) {
                Tools.showMessageDialog(ex1.getMessage());
            }
            return false;
        }
    }

    /**
     * Retrieves the count of records in a table.
     *
     * @param tablename the name of the table
     * @return the count of records
     */
    public static int getCount(String tablename) {
        try {
            setConnection();
            Statement st = con.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
            String auto = "SELECT COUNT(*) AS count FROM " + tablename;
            st.executeQuery(auto);
            String num = "";
            while (st.getResultSet().next()) {
                num = st.getResultSet().getString("count");
            }
            con.close();
            if (num == null || "".equals(num)) {
                return 0;
            } else {
                return Integer.parseInt(num);
                //((Integer)num + 1)
            }
        } catch (SQLException ex) {
            Tools.showMessageDialog(ex.getMessage());
            try {
                con.close();
            } catch (SQLException ex1) {
                Tools.showMessageDialog(ex1.getMessage());
            }
            return 0;
        }
    }

    /**
     * Generates an auto number for a specific table and column.
     *
     * @param tablename the name of the table
     * @param colmnname the name of the column
     * @return the generated auto number
     */
    public static String getAutoNumber(String tablename, String colmnname) {
        try {
            setConnection();
            Statement st = con.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
            String auto = "SELECT MAX(" + colmnname + ")+1 AS Autonum FROM " + tablename;
            st.executeQuery(auto);
            String num = "";
            while (st.getResultSet().next()) {
                num = st.getResultSet().getString("Autonum");
            }
            con.close();
            if (num == null || "".equals(num)) {
                return "1";
            } else {
                return num;
                //((Integer)num + 1)
            }
        } catch (SQLException ex) {
            Tools.showMessageDialog(ex.getMessage());
            try {
                con.close();
            } catch (SQLException ex1) {
                Tools.showMessageDialog(ex1.getMessage());
            }
            return "0";
        }
    }

    /**
     * Retrieves data from the database based on the provided SQL statement.
     *
     * @param statement the SQL statement to execute
     * @return a Table object containing the retrieved data
     */
    public static Table getTableData(String statement) {
        try {
            setConnection();
            Statement st = con.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
            ResultSet rs = st.executeQuery(statement);
            ResultSetMetaData rsmd = rs.getMetaData();
            int c = rsmd.getColumnCount();
            Table t = new Table(c);
            while (rs.next()) {
                Object row[] = new Object[c];
                for (int i = 0; i < c; i++) {
                    row[i] = rs.getString(i + 1);
                }
                t.addNewRow(row);
            }
            con.close();
            return t;
        } catch (SQLException ex) {
            Tools.showMessageDialog(ex.getMessage());
            return new Table(0);
        }
    }

    /**
     * Checks if a value exists in a table column.
     *
     * @param tablename the name of the table
     * @param primarykeycolmnname the name of the primary key column
     * @param value the value to check
     * @return true if the value exists, false otherwise
     */
    public static boolean isExist(String tablename, String primarykeycolmnname, String value) {
        try {
            setConnection();
            Statement st = con.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
            String query = "SELECT " + primarykeycolmnname + " FROM " + tablename;
            ResultSet rs = st.executeQuery(query);
            rs.last();
            int c = rs.getRow();
            rs.beforeFirst();
            String values[] = new String[c];
            int i = 0;
            while (rs.next()) {
                values[i] = rs.getString(1);
                i++;
            }
            for (int j = 0; j < c; j++) {
                if (values[j].equals(value)) {
                    con.close();
                    return true;
                }
            }
            con.close();
            return false;
        } catch (SQLException ex) {
            Tools.showMessageDialog(ex.getMessage());
        }
        return true;
    }

    /**
     * Fills a custom combo box with values from a table column.
     *
     * @param tablename the name of the table
     * @param colmnname the name of the column
     * @param ccb the custom combo box to fill
     */
    public static void fillComboBox(String tablename, String colmnname, CustomComboBox ccb) {
        try {
            setConnection();
            Statement st = con.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
            String query = "SELECT " + colmnname + " FROM " + tablename;
            ResultSet rs = st.executeQuery(query);
            rs.last();
            int c = rs.getRow();
            rs.beforeFirst();
            String values[] = new String[c];
            int i = 0;
            while (rs.next()) {
                values[i] = rs.getString(1);
                i++;
            }
            con.close();
            ccb.setModel(new DefaultComboBoxModel(values));
        } catch (SQLException ex) {
            Tools.showMessageDialog(ex.getMessage());
        }
    }

    /**
     * Fills a table with data based on the provided table name or SQL
     * statement.
     *
     * @param Tname_Statement the table name or SQL statement
     * @param table the table to fill
     */
    public static void fillTable(String Tname_Statement, JTable table) {
        try {
            setConnection();
            Statement st = con.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
            ResultSet rs;
            String query;
            if (Tname_Statement.length() > 8) {
                String isQ = Tname_Statement.substring(0, 7).toLowerCase();
                if (isQ.equals("select ")) {
                    query = Tname_Statement;
                } else {
                    query = "SELECT * FROM " + Tname_Statement;
                }
            } else {
                query = "SELECT * FROM " + Tname_Statement;
            }
            rs = st.executeQuery(query);
            ResultSetMetaData tsmd = rs.getMetaData();
            int c = tsmd.getColumnCount();
            DefaultTableModel dtm = (DefaultTableModel) table.getModel();
            Vector row = new Vector();
            dtm.setRowCount(0);
            while (rs.next()) {
                row = new Vector(c);
                for (int i = 1; i <= c; i++) {
                    row.add(rs.getString(i));
                }
                dtm.addRow(row);
            }
            if (table.getColumnCount() != c) {
                Tools.showMessageDialog("Tables columns count ( " + table.getColumnCount() + " ) not equals to Query columns count ( " + c + " )");
            }
            con.close();
        } catch (SQLException ex) {
            Tools.showMessageDialog(ex.getMessage());
        }
    }

    /**
     * Retrieves a value from a specific column in a table based on the provided
     * primary key value.
     *
     * @param tablename the name of the table
     * @param primarykey the name of the primary key column
     * @param value the primary key value
     * @return the retrieved value
     */
    public static String getValue(String tablename, String primarykey, String value) {
        try {
            setConnection();
            Statement st = con.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE, ResultSet.CONCUR_READ_ONLY);
            ResultSet rs = st.executeQuery("SELECT * FROM " + tablename + " WHERE " + primarykey + " = '" + value + "'");
            rs.next();
            String result = rs.getString(3);
            con.close();
            return result;
        } catch (SQLException ex) {
            Tools.showMessageDialog(ex + "");
        }
        return null;
    }

    /**
     * Initializes the database by executing a series of SQL statements.
     *
     * @return true if the database is successfully initialized, false
     * otherwise.
     */
    public static boolean initializeDatabase() {
        // SQL statements to create tables, views, and other database objects
        String[] sqlStatements = {
            "CREATE DATABASE IF NOT EXISTS Company DEFAULT CHARACTER SET UTF8 DEFAULT COLLATE UTF8_GENERAL_CI",
            "USE Company",
            "CREATE TABLE IF NOT EXISTS Department (DepartmentNumber INT PRIMARY KEY, DepartmentName VARCHAR(150) UNIQUE, Location VARCHAR(100))",
            "CREATE TABLE IF NOT EXISTS Project (ProjectNumber INT PRIMARY KEY, ProjectName VARCHAR(150) UNIQUE, Location VARCHAR(100), DepartmentNumber INT, FOREIGN KEY (DepartmentNumber) REFERENCES Department(DepartmentNumber))",
            "CREATE TABLE IF NOT EXISTS Employee (EmployeeNumber INT PRIMARY KEY, EmployeeName VARCHAR(150) UNIQUE, Address VARCHAR(200), Salary DOUBLE, HiringDate DATE, BirthDate DATE, DepartmentNumber INT, FOREIGN KEY (DepartmentNumber) REFERENCES Department(DepartmentNumber))",
            "CREATE TABLE IF NOT EXISTS EmployeePhones (Phone VARCHAR(20), EmployeeNumber INT, PRIMARY KEY (Phone, EmployeeNumber), FOREIGN KEY (EmployeeNumber) REFERENCES Employee(EmployeeNumber))",
            "CREATE TABLE IF NOT EXISTS WorkOn (EmployeeNumber INT, ProjectNumber INT, PRIMARY KEY (ProjectNumber, EmployeeNumber), FOREIGN KEY (ProjectNumber) REFERENCES Project(ProjectNumber), FOREIGN KEY (EmployeeNumber) REFERENCES Employee(EmployeeNumber))",
            "CREATE TABLE IF NOT EXISTS Users (Username VARCHAR(150) PRIMARY KEY, Password VARCHAR(150), IsBanned varchar(150))",
            "INSERT IGNORE INTO Users VALUES ('admin', 'admin', 'No')",
            "CREATE VIEW IF NOT EXISTS EmployeeData AS SELECT EmployeeNumber AS 'Number', EmployeeName AS 'Name', Address, Salary, HiringDate AS 'HiringDate', BirthDate AS 'BirthDate', Employee.DepartmentNumber, DepartmentName AS 'Department' FROM Department, Employee WHERE Employee.DepartmentNumber = Department.DepartmentNumber",
            "CREATE VIEW IF NOT EXISTS DepartmentData AS SELECT DepartmentNumber, DepartmentName AS 'Department', Location FROM Department",
            "CREATE VIEW IF NOT EXISTS ProjectData AS SELECT ProjectNumber, ProjectName AS 'Project', Project.DepartmentNumber AS 'DepartmentNumber', Department.DepartmentName AS 'Department', Project.Location FROM Project, Department WHERE Project.DepartmentNumber = Department.DepartmentNumber",
            "CREATE VIEW IF NOT EXISTS EmployeePhonesData AS SELECT Employee.EmployeeNumber AS 'EmployeeNumber', Employee.EmployeeName AS 'EmployeeName', EmployeePhones.Phone FROM Employee, EmployeePhones WHERE Employee.EmployeeNumber = EmployeePhones.EmployeeNumber",
            "CREATE VIEW IF NOT EXISTS WorkOnData AS SELECT WorkOn.EmployeeNumber AS 'EmployeeNumber', Employee.EmployeeName AS 'EmployeeName', WorkOn.ProjectNumber AS 'ProjectNumber', Project.ProjectName AS 'ProjectName' FROM WorkOn, Employee, Project WHERE WorkOn.EmployeeNumber = Employee.EmployeeNumber AND WorkOn.ProjectNumber = Project.ProjectNumber"
        };

        try {
            // Establish the database connection
            setConnection();

            // Create a statement for executing SQL queries
            Statement statement = con.createStatement();

            // Execute each SQL statement
            for (String sql : sqlStatements) {
                statement.executeUpdate(sql);
            }

            // Close the database connection
            con.close();
            DBName = "Company?useUnicode=true#&characterEncoding=UTF-8";
            /* Allow Arabic in the database */
            return true;

        } catch (SQLException ex) {
            // Handle any SQL errors that occur during initialization
            Tools.showMessageDialog(ex.getMessage());

            try {
                // Close the database connection in case of an error
                con.close();
            } catch (SQLException ex1) {
                Tools.showMessageDialog(ex1.getMessage());
            }
            return false;
        }
    }

}
