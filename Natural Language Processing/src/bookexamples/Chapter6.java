package bookexamples;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Properties;

import opennlp.tools.doccat.DoccatModel;
import opennlp.tools.doccat.DocumentCategorizerME;
import opennlp.tools.doccat.DocumentSample;
import opennlp.tools.doccat.DocumentSampleStream;
import opennlp.tools.util.ObjectStream;
import opennlp.tools.util.PlainTextByLineStream;

import com.aliasi.classify.Classification;
import com.aliasi.classify.Classified;
import com.aliasi.classify.DynamicLMClassifier;
import com.aliasi.lm.NGramProcessLM;
import com.aliasi.util.Files;

import edu.stanford.nlp.classify.Classifier;
import edu.stanford.nlp.classify.ColumnDataClassifier;
import edu.stanford.nlp.ling.CoreAnnotations;
import edu.stanford.nlp.ling.Datum;
import edu.stanford.nlp.neural.rnn.RNNCoreAnnotations;
import edu.stanford.nlp.objectbank.ObjectBank;
import edu.stanford.nlp.pipeline.Annotation;
import edu.stanford.nlp.pipeline.StanfordCoreNLP;
import edu.stanford.nlp.sentiment.SentimentCoreAnnotations.SentimentAnnotatedTree;
import edu.stanford.nlp.trees.Tree;
import edu.stanford.nlp.util.CoreMap;

/**
 * Chapter 6 examples and tests.
 * 
 * @author Jacob Malter learning from Natural Language Processing with Java
 *
 */
public class Chapter6 {

	@SuppressWarnings("deprecation")
	public static void opennlpClassification() {
		DoccatModel model = null;
		try (InputStream dataIn = new FileInputStream("en-animal.train");
				OutputStream dataOut = new FileOutputStream("en-animal.model");) {
			ObjectStream<String> lineStream = new PlainTextByLineStream(dataIn,
					"UTF-8");
			ObjectStream<DocumentSample> sampleStream = new DocumentSampleStream(
					lineStream);
			model = DocumentCategorizerME.train("en", sampleStream);
			OutputStream modelOut = null;
			modelOut = new BufferedOutputStream(dataOut);
			model.serialize(modelOut);
			System.out.println();
		} catch (IOException ex) {
			System.out.println("IO");
		}
	}

	public static void opennlpDocumentCategorizer() {
		final String toto = "Toto belongs to Dorothy Gale, the heroine of "
				+ "the first and many subsequent books. In the first "
				+ "book, he never spoke, although other animals, native "
				+ "to Oz, did. In subsequent books, other animals "
				+ "gained the ability to speak upon reaching Oz or "
				+ "similar lands, but Toto remained speechless.";
		try (InputStream modelIn = new FileInputStream(new File(
				"en-animal.model"));) {
			DoccatModel model = new DoccatModel(modelIn);
			DocumentCategorizerME categorizer = new DocumentCategorizerME(model);
			double outcomes[] = categorizer.categorize(toto);
			for (int i = 0; i < categorizer.getNumberOfCategories(); i++) {
				String category = categorizer.getCategory(i);
				System.out.println(category + " - " + outcomes[i]);
			}
		} catch (IOException ex) {
			System.out.println("IO");
		}
	}

	public static void stanfordClassification() {
		ColumnDataClassifier cdc = new ColumnDataClassifier("box.prop");
		Classifier<String, String> classifier = cdc.makeClassifier(cdc
				.readTrainingExamples("box.train"));
		for (String line : ObjectBank.getLineIterator("box.test", "utf-8")) {
			Datum<String, String> datum = cdc.makeDatumFromLine(line);
			System.out.println("Datum: {" + line + "]\tPredicted Category: "
					+ classifier.classOf(datum));
		}
	}

	@SuppressWarnings("unused")
	public static void stanfordPipeline() {
		String review = "An overly sentimental film with a somewhat "
				+ "problematic message, but its sweetness and charm "
				+ "are occasionally enough to approximate true depth "
				+ "and grace. ";
		String sam = "Sam was an odd sort of fellow. Not prone "
				+ "to angry and not prone to merriment. Overall, "
				+ "an odd fellow.";
		String mary = "Mary thought that custard pie was the "
				+ "best pie in the world. However, she loathed "
				+ "chocolate pie.";
		Properties props = new Properties();
		props.put("annotators", "tokenize, ssplit, parse, sentiment");
		StanfordCoreNLP pipeline = new StanfordCoreNLP(props);
		Annotation annotation = new Annotation(review);
		pipeline.annotate(annotation);
		String sentimentText[] = { "Very Negative", "Negative", "Neutral",
				"Positive", "Very Positive" };
		for (CoreMap sentence : annotation
				.get(CoreAnnotations.SentencesAnnotation.class)) {
			Tree tree = sentence.get(SentimentAnnotatedTree.class);
			int score = RNNCoreAnnotations.getPredictedClass(tree);
			System.out.println(sentimentText[score]);
		}
	}

	public static void lingpipeClassify() {
		String categories[] = { "soc.religion.christian", "talk.religion.misc",
				"alt.atheism", "misc.forsale" };
		int nGramsize = 6;
		DynamicLMClassifier<NGramProcessLM> classifier = DynamicLMClassifier
				.createNGramProcess(categories, nGramsize);
		String directory = "add file directory/demos";
		File trainingDirectory = new File(directory
				+ "/data/fourNewsGroups/4news-train");
		for (int i = 0; i < categories.length; ++i) {
			File classDir = new File(trainingDirectory, categories[i]);
			String trainingFiles[] = classDir.list();
			for (int j = 0; trainingFiles != null && j < trainingFiles.length; ++j) {
				try {
					File file = new File(classDir, trainingFiles[j]);
					String text = Files.readFromFile(file, "ISO-8859-1");
					Classification classification = new Classification(
							categories[i]);
					Classified<CharSequence> classified = new Classified<>(
							text, classification);
					classifier.handle(classified);
				} catch (IOException ex) {
					System.out.println("IO " + j);
				}
			}
		}
	}

	/**
	 * Runs chapter 6 examples.
	 * 
	 * @param args
	 *            command line arguments not used in this program
	 */
	public static void main(String[] args) {
		opennlpClassification();
		opennlpDocumentCategorizer();
		// stanfordClassification();
		// stanfordPipeline();
		lingpipeClassify();
	}

}