package bookexamples;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Reader;
import java.io.StringReader;
import java.text.BreakIterator;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.aliasi.chunk.Chunk;
import com.aliasi.chunk.Chunking;
import com.aliasi.sentences.IndoEuropeanSentenceModel;
import com.aliasi.sentences.MedlineSentenceModel;
import com.aliasi.sentences.SentenceChunker;
import com.aliasi.tokenizer.IndoEuropeanTokenizerFactory;
import com.aliasi.tokenizer.Tokenizer;
import com.aliasi.tokenizer.TokenizerFactory;

import edu.stanford.nlp.ling.CoreLabel;
import edu.stanford.nlp.ling.HasWord;
import edu.stanford.nlp.pipeline.Annotation;
import edu.stanford.nlp.pipeline.StanfordCoreNLP;
import edu.stanford.nlp.process.CoreLabelTokenFactory;
import edu.stanford.nlp.process.DocumentPreprocessor;
import edu.stanford.nlp.process.PTBTokenizer;
import edu.stanford.nlp.process.WordToSentenceProcessor;
import opennlp.tools.sentdetect.SentenceDetectorME;
import opennlp.tools.sentdetect.SentenceModel;
import opennlp.tools.sentdetect.SentenceSample;
import opennlp.tools.sentdetect.SentenceSampleStream;
import opennlp.tools.util.ObjectStream;
import opennlp.tools.util.PlainTextByLineStream;
import opennlp.tools.util.Span;
import opennlp.tools.util.TrainingParameters;

/**
 * Chapter 3 examples and tests.
 * 
 * @author Jacob Malter learning from Natural Language Processing with Java
 *
 */
public class Chapter3 {

	private static String paragraph = "When determining the end of sentences "
			+ "we need to consider several factors. Sentences may end with "
			+ "exclamation marks! Or possibly questions marks? Within "
			+ "sentences we may find numbers like 3.14159 , abbreviations "
			+ "such as found in Mr. Smith, and possibly ellipses either "
			+ "within a sentence … , or at the end of a sentence…";

	public static void javaSplitSBD() {
		String simple = "[.?!]";
		String[] splitString = paragraph.split(simple);
		for (String string : splitString)
			System.out.println(string);
	}

	public static void javaPatternSBD() {
		Pattern sentencePattern = Pattern
				.compile(
						"# Match a sentence ending in punctuation or EOS.\n"
								+ "[^.!?\\s]    # First char is non-punct, non-ws\n"
								+ "[^.!?]*      # Greedily consume up to punctuation.\n"
								+ "(?:          # Group for unrolling the loop.\n"
								+ "  [.!?]      # (special) inner punctuation ok if\n"
								+ "  (?!['\"]?\\s|$)  # not followed by ws or EOS.\n"
								+ "  [^.!?]*    # Greedily consume up to punctuation.\n"
								+ ")*           # Zero or more (special normal*)\n"
								+ "[.!?]?       # Optional ending punctuation.\n"
								+ "['\"]?       # Optional closing quote.\n"
								+ "(?=\\s|$)", Pattern.MULTILINE
								| Pattern.COMMENTS);
		Matcher matcher = sentencePattern.matcher(paragraph);
		while (matcher.find())
			System.out.println(matcher.group());
	}

	public static void javaBreakIterator() {
		BreakIterator sentenceIterator = BreakIterator.getSentenceInstance();
		sentenceIterator.setText(paragraph);
		int boundary = sentenceIterator.first();
		while (boundary != BreakIterator.DONE) {
			int begin = boundary;
			System.out.print(boundary + "-");
			boundary = sentenceIterator.next();
			int end = boundary;
			if (end == BreakIterator.DONE)
				break;
			System.out.println(boundary + " ["
					+ paragraph.substring(begin, end) + "]");
		}
	}

	public static void openNLPSBD() {
		try (InputStream is = new FileInputStream(new File(
				"add file directory", "en-sent.bin"))) {
			SentenceModel model = new SentenceModel(is);
			SentenceDetectorME detector = new SentenceDetectorME(model);
			String sentences[] = detector.sentDetect(paragraph);
			for (String sentence : sentences)
				System.out.println(sentence);
			System.out.println();
			double probabilities[] = detector.getSentenceProbabilities();
			for (double probability : probabilities)
				System.out.println(probability);
			System.out.println();
			Span spans[] = detector.sentPosDetect(paragraph);
			for (Span span : spans)
				System.out.println(span + "["
						+ paragraph.substring(span.getStart(), span.getEnd())
						+ "]");
			System.out.println();
			paragraph = " This sentence starts with spaces and ends with "
					+ "spaces . This sentences has no spaces between the next "
					+ "one.This is the next one.";
			sentences = detector.sentDetect(paragraph);
			for (String sentence : sentences)
				System.out.println(sentence);
		} catch (FileNotFoundException ex) {
			System.out.println("FileNotFound");
		} catch (IOException ex) {
			System.out.println("IO");
		}
		paragraph = "When determining the end of sentences "
				+ "we need to consider several factors. Sentences may end with "
				+ "exclamation marks! Or possibly questions marks? Within "
				+ "sentences we may find numbers like 3.14159 , abbreviations "
				+ "such as found in Mr. Smith, and possibly ellipses either "
				+ "within a sentence … , or at the end of a sentence…";
	}

	public static void standfordSBD() {
		PTBTokenizer<CoreLabel> ptb = new PTBTokenizer<>(new StringReader(
				paragraph), new CoreLabelTokenFactory(), null);
		WordToSentenceProcessor<CoreLabel> wtsp = new WordToSentenceProcessor<>();
		List<List<CoreLabel>> sents = wtsp.process(ptb.tokenize());
		for (List<CoreLabel> sent : sents) {
			for (CoreLabel element : sent) {
				System.out.print(element + " ");
			}
			System.out.println();
		}
		System.out.println();
		for (List<CoreLabel> sent : sents) {
			for (CoreLabel element : sent) {
				System.out.print(element.endPosition() + " ");
			}
			System.out.println();
		}
		System.out.println();
		for (List<CoreLabel> sent : sents)
			System.out.println(sent.get(0) + " " + sent.get(0).beginPosition());
		System.out.println();
		for (List<CoreLabel> sent : sents) {
			int size = sent.size();
			System.out.println(sent.get(size - 1) + " "
					+ sent.get(size - 1).endPosition());
		}
	}

	public static void standfordDocument() {
		Reader reader = new StringReader(paragraph);
		DocumentPreprocessor dp = new DocumentPreprocessor(reader);
		for (List<HasWord> sentence : dp)
			System.out.println(sentence);
		System.out.println();
		try {
			reader = new FileReader("add file directory");
			dp = new DocumentPreprocessor(reader,
					DocumentPreprocessor.DocType.XML);
			dp.setElementDelimiter("sentence");
			for (List<HasWord> sentence : dp)
				System.out.println(sentence);
		} catch (FileNotFoundException ex) {
			System.out.println("FileNotFound");
		}
	}

	public static void standfordCoreNLP() {
		Properties properties = new Properties();
		properties.put("annotators", "tokenize, ssplit");
		StanfordCoreNLP pipeline = new StanfordCoreNLP(properties);
		Annotation annotation = new Annotation(paragraph);
		pipeline.annotate(annotation);
	}

	public static void lingpipeSBD() {
		// add period at end
		paragraph = "When determining the end of sentences "
				+ "we need to consider several factors. Sentences may end with "
				+ "exclamation marks! Or possibly questions marks? Within "
				+ "sentences we may find numbers like 3.14159 , abbreviations "
				+ "such as found in Mr. Smith, and possibly ellipses either "
				+ "within a sentence … , or at the end of a sentence….";
		TokenizerFactory TOKENIZER_FACTORY = IndoEuropeanTokenizerFactory.INSTANCE;
		IndoEuropeanSentenceModel sentenceModel = new IndoEuropeanSentenceModel();
		List<String> tokenList = new ArrayList<>();
		List<String> whiteList = new ArrayList<>();
		Tokenizer tokenizer = TOKENIZER_FACTORY.tokenizer(
				paragraph.toCharArray(), 0, paragraph.length());
		tokenizer.tokenize(tokenList, whiteList);
		String tokens[] = new String[tokenList.size()];
		String whites[] = new String[whiteList.size()];
		tokenList.toArray(tokens);
		whiteList.toArray(whites);
		int sentenceBoundaries[] = sentenceModel
				.boundaryIndices(tokens, whites);
		int start = 0;
		for (int boundary : sentenceBoundaries) {
			while (start <= boundary) {
				System.out.print(tokenList.get(start)
						+ whiteList.get(start + 1));
				start++;
			}
			System.out.println();
		}
		// remove period at end
		paragraph = "When determining the end of sentences "
				+ "we need to consider several factors. Sentences may end with "
				+ "exclamation marks! Or possibly questions marks? Within "
				+ "sentences we may find numbers like 3.14159 , abbreviations "
				+ "such as found in Mr. Smith, and possibly ellipses either "
				+ "within a sentence … , or at the end of a sentence…";
	}

	public static void lingpipeChunker() {
		paragraph = "When determining the end of sentences "
				+ "we need to consider several factors. Sentences may end with "
				+ "exclamation marks! Or possibly questions marks? Within "
				+ "sentences we may find numbers like 3.14159 , abbreviations "
				+ "such as found in Mr. Smith, and possibly ellipses either "
				+ "within a sentence … , or at the end of a sentence….";
		TokenizerFactory tokenizerfactory = IndoEuropeanTokenizerFactory.INSTANCE;
		IndoEuropeanSentenceModel sentenceModel = new IndoEuropeanSentenceModel();
		SentenceChunker sentenceChunker = new SentenceChunker(tokenizerfactory,
				sentenceModel);
		Chunking chunking = sentenceChunker.chunk(paragraph.toCharArray(), 0,
				paragraph.length());
		Set<Chunk> sentences = chunking.chunkSet();
		String slice = chunking.charSequence().toString();
		for (Chunk sentence : sentences)
			System.out.println("["
					+ slice.substring(sentence.start(), sentence.end()) + "]");
		paragraph = "When determining the end of sentences "
				+ "we need to consider several factors. Sentences may end with "
				+ "exclamation marks! Or possibly questions marks? Within "
				+ "sentences we may find numbers like 3.14159 , abbreviations "
				+ "such as found in Mr. Smith, and possibly ellipses either "
				+ "within a sentence … , or at the end of a sentence…";
	}

	public static void lingpipeMedline() {
		paragraph = "HepG2 cells were obtained from the American Type Culture "
				+ "Collection (Rockville, MD, USA) and were used only until "
				+ "passage 30. They were routinely grown at 37°C in Dulbecco's "
				+ "modified Eagle's medium (DMEM) containing 10 % fetal bovine "
				+ "serum (FBS), 2 mM glutamine, 1 mM sodium pyruvate, and 25 "
				+ "mM glucose (Invitrogen, Carlsbad, CA, USA) in a humidified "
				+ "atmosphere containing 5% CO2. For precursor and 13C-sugar "
				+ "experiments, tissue culture treated polystyrene 35 mm "
				+ "dishes (Corning Inc, Lowell, MA, USA) were seeded with 2 "
				+ "× 106 cells and grown to confluency in DMEM.";
		TokenizerFactory tokenizerfactory = IndoEuropeanTokenizerFactory.INSTANCE;
		MedlineSentenceModel sentenceModel = new MedlineSentenceModel();
		SentenceChunker sentenceChunker = new SentenceChunker(tokenizerfactory,
				sentenceModel);
		Chunking chunking = sentenceChunker.chunk(paragraph.toCharArray(), 0,
				paragraph.length());
		Set<Chunk> sentences = chunking.chunkSet();
		String slice = chunking.charSequence().toString();
		for (Chunk sentence : sentences)
			System.out.println("["
					+ slice.substring(sentence.start(), sentence.end()) + "]");
		paragraph = "When determining the end of sentences "
				+ "we need to consider several factors. Sentences may end with "
				+ "exclamation marks! Or possibly questions marks? Within "
				+ "sentences we may find numbers like 3.14159 , abbreviations "
				+ "such as found in Mr. Smith, and possibly ellipses either "
				+ "within a sentence … , or at the end of a sentence…";
	}

	public static void trainingSentenceModel() {
		try {
			@SuppressWarnings("deprecation")
			ObjectStream<String> lineStream = new PlainTextByLineStream(
					new FileReader("sentence.train"));
			ObjectStream<SentenceSample> sampleStream = new SentenceSampleStream(
					lineStream);
			@SuppressWarnings("deprecation")
			SentenceModel model = SentenceDetectorME.train("en", sampleStream,
					true, null, TrainingParameters.defaultParams());
			OutputStream modelStream = new BufferedOutputStream(
					new FileOutputStream("modelFile"));
			model.serialize(modelStream);
		} catch (FileNotFoundException ex) {
			System.out.println("FileNotFound");
		} catch (IOException ex) {
			System.out.println("IO");
		}
	}

	/**
	 * Runs chapter 3 examples.
	 * 
	 * @param args
	 *            command line arguments not used in this program
	 */
	public static void main(String[] args) {
		javaSplitSBD();
		System.out.println();
		javaPatternSBD();
		System.out.println();
		javaBreakIterator();
		System.out.println();
		openNLPSBD();
		System.out.println();
		standfordSBD();
		System.out.println();
		standfordDocument();
		System.out.println();
		standfordCoreNLP();
		System.out.println();
		lingpipeSBD();
		System.out.println();
		lingpipeChunker();
		System.out.println();
		lingpipeMedline();
		System.out.println();
		trainingSentenceModel();
	}

}