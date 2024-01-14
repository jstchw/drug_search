import {Container, Row, Col} from 'react-bootstrap'
import Header from '../../components/Header/Header'
import SearchBox from '../../components/SearchBox/SearchBox'
import './Search.css'


const Search = () => {
    // const [eventResults, setEventResults] = React.useState({})
    // const [eventsOverTime, setEventsOverTime] = React.useState([])
    //
    //
    // const [isLoading, setIsLoading] = React.useState(false)
    // const [searchError, setSearchError] = React.useState(false)
    //
    // const [currentSearchTerm, setCurrentSearchTerm] = React.useState('')
    // const [mainSearchTerm, setMainSearchTerm] = React.useState('')
    // const [additionalSearchTerm, setAdditionalSearchTerm] = React.useState('')
    //
    // const [showAdditionalSearch, setShowAdditionalSearch] = React.useState(false)
    //
    // const [selectedSearchOptionIndex, setSelectedSearchOptionIndex] = React.useState({
    //     searchBy: searchTypes.findIndex(({ value }) => value === Cookies.get('searchBy')) ?? searchTypes[0]?.index,
    //     sex: searchSex.findIndex(({ value }) => value === Cookies.get('lastSelectedSex')) ?? searchSex[0]?.index,
    //     country: searchCountry.findIndex(({ value }) => value === Cookies.get('lastSelectedCountry') ?? 'US'),
    // })


    // This has to be changed to a whole object instead values
    // const [searchOptions, setSearchOptions] = React.useState({
    //     searchBy: searchTypes[selectedSearchOptionIndex.searchBy].value,
    //     sex: searchSex[selectedSearchOptionIndex.sex],
    //     age: [
    //         {
    //             value: Cookies.get('lastSelectedAgeStart') ?? searchAgeRange[0],
    //             index: 0,
    //             label: 'Minimum Age',
    //             type: 'age'
    //         },
    //         {
    //             value: Cookies.get('lastSelectedAgeEnd') ?? searchAgeRange[1],
    //             index: 1,
    //             label: 'Maximum Age',
    //             type: 'age'
    //         },
    //     ],
    //     country: searchCountry[selectedSearchOptionIndex.country],
    // })

    // const [searchHistory, setSearchHistory] = React.useState(JSON.parse(Cookies.get('searchHistory') || '[]'))


    // useEffect(() => {
    //     if(!showAdditionalSearch) {
    //         setAdditionalSearchTerm('')
    //     }
    // }, [showAdditionalSearch])

    // Manipulate the cookies and display of drug search terms
    // useEffect(() => {
    //     if (currentSearchTerm) {
    //         const capitalizedTerms = currentSearchTerm.map((term: string) => term.charAt(0).toUpperCase() + term.slice(1));
    //         document.title = capitalizedTerms.join(' & ') + ' - DrugSearch';
    //         updateSearchHistory(capitalizedTerms, searchHistory, setSearchHistory, searchOptions)
    //     }
    // }, [currentSearchTerm]);


    // const handleSearch = async (searchOptions, singleSearchTerm) => {
    //
    //     let events
    //
    //     setEventResults(null)
    //     setSearchError(false)
    //
    //
    //     // This is a very, very bad implementation but it works
    //     // 100% needs to be refactored
    //     // When a search is submitted from a cookie button - the search term is an array from the beginning
    //     // When a search is submitted from the search box - the search term needs to be formed into an array and checked if other search term is present
    //     if(singleSearchTerm) {
    //         setIsLoading(true)
    //
    //         let searchTerm
    //         if(typeof singleSearchTerm === 'string') {
    //             searchTerm = additionalSearchTerm ? [mainSearchTerm, additionalSearchTerm] : [singleSearchTerm]
    //         } else {
    //             searchTerm = singleSearchTerm
    //         }
    //
    //         if (mainSearchTerm === additionalSearchTerm) {
    //             setSearchError(true)
    //             setIsLoading(false)
    //             return
    //         }
    //
    //         console.log(searchOptions)
    //         // Search for an ADE and get drugs
    //         if (searchOptions.searchBy.value === searchTypes[2]?.value) {
    //             events = await getDrugsFromEvents(searchTerm, searchOptions)
    //         } else {
    //             events = await getEventsFromDrugs(searchTerm, searchOptions)
    //         }
    //
    //         const eventsOverTime = await getEventsOverTime(searchTerm, searchOptions)
    //
    //         setIsLoading(false)
    //
    //         if(events.result) {
    //             setEventResults(events)
    //             setCurrentSearchTerm(searchTerm)
    //         } else {
    //             setSearchError(true)
    //         }
    //
    //         if(eventsOverTime.result) {
    //             setEventsOverTime(eventsOverTime)
    //         } else {
    //             setSearchError(true)
    //         }
    //     }
    // }

    return (
        <div className="min-vh-100 d-flex flex-column">
            <Container fluid className={'d-flex flex-grow-1 align-items-center'}>
                <Row className="w-100 mx-auto">
                    <Col xs="auto" className="text-center mb-5 mt-3 mx-auto">
                        <Row className={'mb-4'}>
                            <Container className={'mb-4 mt-3 w-auto'}>
                                <Header />
                            </Container>
                            <SearchBox/>
                        </Row>
                    </Col>
                </Row>
            </Container>
        </div>
    );
}

export default Search;