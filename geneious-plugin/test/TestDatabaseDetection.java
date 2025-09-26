import com.jhuapl.taxtriage.geneious.database.DatabaseManager;
import com.jhuapl.taxtriage.geneious.database.DatabaseManager.DatabaseType;

/**
 * Simple test to verify database detection logic
 */
public class TestDatabaseDetection {

    public static void main(String[] args) {
        System.out.println("Testing TaxTriage Database Detection");
        System.out.println("====================================\n");

        DatabaseManager dbManager = DatabaseManager.getInstance();

        // Test viral database detection
        System.out.println("Testing VIRAL database:");
        System.out.println("-----------------------");

        // First test without allowing download
        String viralPath = dbManager.getDatabasePath(DatabaseType.VIRAL, false);
        if (viralPath != null) {
            System.out.println("✓ Cached viral database found at: " + viralPath);
            System.out.println("  Expected behavior: Use cached database with --download_db false");
        } else {
            System.out.println("✗ No cached viral database found");
            System.out.println("  Expected behavior: Should download with --download_db true");
        }

        // Test with allowing download
        String viralPathWithDownload = dbManager.getDatabasePath(DatabaseType.VIRAL, true);
        if (viralPathWithDownload != null && !"DOWNLOAD_REQUIRED".equals(viralPathWithDownload)) {
            System.out.println("✓ Using existing database: " + viralPathWithDownload);
        } else if ("DOWNLOAD_REQUIRED".equals(viralPathWithDownload)) {
            System.out.println("⚠ Download required for viral database");
        }

        System.out.println("\nTesting STANDARD database:");
        System.out.println("-------------------------");

        // Test standard database
        String standardPath = dbManager.getDatabasePath(DatabaseType.STANDARD, false);
        if (standardPath != null) {
            System.out.println("✓ Cached standard database found at: " + standardPath);
            System.out.println("  Expected behavior: Use cached database with --download_db false");
        } else {
            System.out.println("✗ No cached standard database found");
            System.out.println("  Expected behavior: Should download with --download_db true");
        }

        System.out.println("\nSummary:");
        System.out.println("--------");
        System.out.println("The database detection logic should:");
        System.out.println("1. Use cached databases when available (--download_db false)");
        System.out.println("2. Download databases when not cached (--download_db true)");
        System.out.println("3. Both hash.k2d AND taxo.k2d must exist for a valid database");
    }
}